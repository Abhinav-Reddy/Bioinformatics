#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <vector>
#include <map>
#include <algorithm>
#include <ctime>
#include <string>

using namespace std;
vector< pair<string, string> > queries;
vector< pair<string, string> > database;
string alphabets;
int mapAlphabets[256];
vector< vector<int> >score;
int neighbors;
int gapPenality;

class AlignmentResult {
	public:
	int score; 
	string ids[2];
	int startPos[2];
	string sequence[2];
};

vector< pair<string, string> > readSequencesFromFile(char* fileName) {
	vector< pair<string, string> > dest;
	ifstream in(fileName);
	string input;
	string sequence = "";
	string id = "";
	while(getline(in, input)) {
		if (input[0] == '>') {
			if (id != "") {
				dest.push_back( make_pair(id, sequence) );
				id = "";
				sequence = "";
				
			}
			
			int i;
			for (i=0;i<input.length() && input[i] != ' '; i++);
			id = input.substr(5, i-5);
		}
		else {
			sequence.append(input); 
		}
	}
	
	dest.push_back( make_pair(id, sequence) );
	return dest;
}


void parseAlphabets(char* fileName) {
	ifstream in(fileName);
	getline(in, alphabets);
	for (int i=0;i<256;i++)
		mapAlphabets[i]=-1;
	for (int i=0;i<alphabets.length();i++) {
		mapAlphabets[alphabets[i]] = i;
		mapAlphabets[alphabets[i]+32] = i;
		mapAlphabets[alphabets[i]-32] = i;
	}
}


void parseScoreMatrix(char* fileName) {
	ifstream in(fileName);
	int tmp;
	for (int i=0;i<alphabets.length();i++) {
		score.push_back(vector<int>(alphabets.size()));
		for (int j=0;j<alphabets.length();j++) {
			in >> tmp;
			score[i][j] = tmp;
		}
	}
	
}


void debugOutput(vector < vector<int> >& dp, AlignmentResult& res) {
	for (int i=0;i<dp.size();i++) {
		for (int j=0;j<dp[0].size();j++)
			cout<<dp[i][j]<<" ";
		cout<<endl;
	}
}

void getAlignmentResult(vector < vector<int> >& dp, vector < vector<int> >& prev,
									int row, int col, string& a, string& b, AlignmentResult& res) {
	res.sequence[0] = "";
	res.sequence[1] = "";
	res.score = dp[row][col];
	while(prev[row][col] != -1) {
		int r = prev[row][col]/(b.size()+1);
		int c = prev[row][col]%(b.size()+1);
		if (r == row-1 && c == col-1) {
			res.sequence[0].push_back(a[row-1]);
			res.sequence[1].push_back(b[col-1]);
		}
		else if (r == row) {
			res.sequence[0].push_back('.');
			res.sequence[1].push_back(b[col-1]);
		}
		else {
			res.sequence[0].push_back(a[row-1]);
			res.sequence[1].push_back('.');
		}
		row = r;
		col = c;
	}
	reverse(res.sequence[0].begin(), res.sequence[0].end());
	reverse(res.sequence[1].begin(), res.sequence[1].end());
	res.startPos[0] = row;
	res.startPos[1] = col;
	//debugOutput(dp, res);
	dp.clear();
	prev.clear();
}


void intializeTables(vector < vector<int> >& dp, vector < vector<int> >& prev, int row, int col) {
	 
}

void getMaxNeighbor(int &p, int& q, int& val, int i, int j, string& a, string& b, vector< vector<int> >& dp) {
	p=i-1;
	q=j-1;
	val = dp[p][q]+score[mapAlphabets[a[p]]][mapAlphabets[b[q]]];
	if (val < dp[i-1][j]+gapPenality) {
		val = dp[i-1][j]+gapPenality;
		p = i-1;
		q = j;
	}
	if (val < dp[i][j-1]+gapPenality) {
		val = dp[i][j-1]+gapPenality;
		p = i;
		q = j-1;
	}
}

void getGlobalAlignment(string& a, string& b, AlignmentResult& res) {
	int row = a.length()+1;
	int col = b.length()+1;
	vector < vector<int> >dp;
	vector < vector<int> > prev;
	dp.push_back(vector<int>(col));
	prev.push_back(vector<int>(col));
	dp[0][0] = 0;
	prev[0][0] = -1;
	for (int i=1;i<row;i++) {
		dp.push_back(vector<int>(col));
		prev.push_back(vector<int>(col));
		dp[i][0] = dp[i-1][0] + gapPenality;
		prev[i][0] = (i-1)*col;
	}
	
	for (int i=1;i<col;i++) {
		dp[0][i] = dp[0][i-1] + gapPenality;
		prev[0][i] = i-1;
	}
	
	
	for (int i=1;i<row;i++) {
		for (int j=1;j<col;j++) {
			int p, q, val;
			getMaxNeighbor(p, q, val, i, j, a, b, dp);
			dp[i][j] = val;
			prev[i][j] = p*col+q;
		}
	}
	getAlignmentResult(dp, prev, row-1, col-1, a, b, res);
}



void getLocalAlignment(string a, string b, AlignmentResult& res) {
	int row = a.size()+1;
	int col = b.size()+1;
	vector < vector<int> >dp;
	vector < vector<int> > prev;
	for (int i=0;i<row;i++) {
		dp.push_back(vector<int>(col));
		prev.push_back(vector<int>(col));
		dp[i][0] = 0;
		prev[i][0] = -1;
	}
	
	for (int i=1;i<col;i++) {
		dp[0][i] = 0;
		prev[0][i] = -1;
	}
	int mx[2];
	mx[0] = mx[1] = 0;
	for (int i=1;i<row;i++) {
		for (int j=1;j<col;j++) {
			int p, q, val;
			getMaxNeighbor(p, q, val, i, j, a, b, dp);
			
			if (val <= 0) {
				dp[i][j] = 0;
				prev[i][j] = -1;
			}
			else {
				dp[i][j] = val;
				prev[i][j] = p*col+q;
			}
			if (dp[mx[0]][mx[1]] < dp[i][j]) {
				mx[0] = i;
				mx[1] = j;
			}
		}
	}
	getAlignmentResult(dp, prev, mx[0], mx[1], a, b, res);
}

void getEndSpaceAlignment(string a, string b, AlignmentResult& res) {
	int row = a.size()+1;
	int col = b.size()+1;
	
	vector < vector<int> >dp;
	vector < vector<int> > prev;
	for (int i=0;i<row;i++) {
		dp.push_back(vector<int>(col));
		prev.push_back(vector<int>(col));
		dp[i][0] = 0;
		prev[i][0] = -1;
	}
	
	for (int i=1;i<col;i++) {
		dp[0][i] = 0;
		prev[0][i] = -1;
	}
	
	for (int i=1;i<row;i++) {
		for (int j=1;j<col;j++) {
			int p, q, val;
			getMaxNeighbor(p, q, val, i, j, a, b, dp);
			
			dp[i][j] = val;
			prev[i][j] = p*col+q;
		}
	}
	int mx[2];
	mx[0] = mx[1] = 0;
	for (int i=1;i<row;i++)
		if (dp[mx[0]][mx[1]] < dp[i][col-1]) {
			mx[0] = i;
			mx[1] = col-1;
		}
		
	for (int i=1;i<col;i++)
		if (dp[mx[0]][mx[1]] < dp[row-1][i]) {
			mx[0] = row-1;
			mx[1] = i;
		}
	getAlignmentResult(dp, prev, mx[0], mx[1], a, b, res);
}

void topKalignments(int alignmentMethod) {
	map<int, vector<AlignmentResult> >mp;
	map<int, vector<int> >execMap;
	for (int i=0;i<queries.size();i++) {
		int start_s=clock();
		for (int j=0;j<database.size();j++) {
			AlignmentResult tmp;
			if (alignmentMethod == 1)
				getGlobalAlignment(queries[i].second, database[j].second, tmp);
			else if (alignmentMethod == 2)
				getLocalAlignment(queries[i].second, database[j].second, tmp);
			else {
				getEndSpaceAlignment(queries[i].second, database[j].second, tmp);
			}
			tmp.ids[1] = database[j].first;
			tmp.ids[0] = queries[i].first;
			//cout<<tmp.score<<endl;
			mp[tmp.score].push_back(tmp);
		}
		start_s = (clock() - start_s)/double(CLOCKS_PER_SEC)*1000;
		execMap[queries[i].second.length()].push_back(start_s);
	}
	int cnt = 0;
	for (map<int, vector<AlignmentResult> >::reverse_iterator  itr=mp.rbegin(); itr!= mp.rend() && cnt < neighbors;itr++) {
		for (int i=0;i<itr->second.size() && cnt < neighbors; i++) {
			cnt++;
			AlignmentResult tmp = itr->second[i];
			cout<<"Score = "<<tmp.score<<endl;
			cout<<tmp.ids[0]<<" "<<tmp.startPos[0]<<" "<<tmp.sequence[0]<<endl;
			cout<<tmp.ids[1]<<" "<<tmp.startPos[1]<<" "<<tmp.sequence[1]<<endl;
		}
	}
}


void cleanSequences() {
	for (int i=0;i<queries.size();i++) {
		string tmpStr = "";
		for (int j=0;j<queries[i].second.length();j++) {
			if (mapAlphabets[queries[i].second[j]] != -1)
				tmpStr.push_back(queries[i].second[j]);
		}
		queries[i].second = tmpStr;
	}
	
	for (int i=0;i<database.size();i++) {
		string tmpStr = "";
		for (int j=0;j<database[i].second.length();j++) {
			if (mapAlphabets[database[i].second[j]] != -1)
				tmpStr.push_back(database[i].second[j]);
		}
		database[i].second = tmpStr;
	}
}



void getGraphData() {
	neighbors = 100;
	string tmpAlignment[4];
	tmpAlignment[1] = "Global";
	tmpAlignment[2] = "Local";
	tmpAlignment[3] = "EndSpace";
	for (int alignmentMethod=1;alignmentMethod<4;alignmentMethod++) {
		map<int, vector<AlignmentResult> >mp;
		ofstream lenFile, timeFile, scoreFile;
  		lenFile.open ((tmpAlignment[alignmentMethod]+"Length.txt").c_str());
  		timeFile.open ((tmpAlignment[alignmentMethod]+"Time.txt").c_str());
		scoreFile.open ((tmpAlignment[alignmentMethod]+"Score.txt").c_str());
		for (int i=0;i<queries.size();i++) {
			
			int start_s=clock();
			for (int j=0;j<database.size();j++) {
				AlignmentResult tmp;
				if (alignmentMethod == 1)
					getGlobalAlignment(queries[i].second, database[j].second, tmp);
				else if (alignmentMethod == 2)
					getLocalAlignment(queries[i].second, database[j].second, tmp);
				else {
					getEndSpaceAlignment(queries[i].second, database[j].second, tmp);
				}
				tmp.ids[1] = database[j].first;
				tmp.ids[0] = queries[i].first;
				//cout<<tmp.score<<endl;
				mp[tmp.score].push_back(tmp);
			}
			start_s = (clock() - start_s)/double(CLOCKS_PER_SEC)*1000000;
			lenFile << queries[i].second.length()<<"\n";
			timeFile << start_s<<"\n";
		}
		
		
		int cnt = 0;
		for (map<int, vector<AlignmentResult> >::reverse_iterator  itr=mp.rbegin(); itr!= mp.rend() && cnt < neighbors;itr++) {
			for (int i=0;i<itr->second.size() && cnt < neighbors; i++) {
				cnt++;
				AlignmentResult tmp = itr->second[i];
				scoreFile<<tmp.score<<"\n";
			}
		}
		lenFile.close();
		timeFile.close();
		scoreFile.close();
	}
}

int main(int argc, char *argv[]) {
	
	if (argc != 8) {
		cout<<"Invalid input"<<endl;
		return 0;
	}
	queries = readSequencesFromFile(argv[2]);
	database = readSequencesFromFile(argv[3]);
	parseAlphabets(argv[4]);
	parseScoreMatrix(argv[5]);
	neighbors = atoi(argv[6]);
	gapPenality = atoi(argv[7]);
	cleanSequences();
	if (atoi(argv[1]) == 4) {
		getGraphData();
	}
	else {
		topKalignments(atoi(argv[1]));
	}
    return 0;
}    
