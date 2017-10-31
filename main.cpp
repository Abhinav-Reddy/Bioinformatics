#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <vector>
#include <map>
#include <ctime>

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
	void operator= (const AlignmentResult& c) {
		startPos[0] = c.startPos[0];
		startPos[1] = c.startPos[1];
		
		sequence[0] = c.sequence[0];
		sequence[1] = c.sequence[1];
		
		score = c.score;
	}
};

void readSequencesFromFile(char* fileName, vector< pair<string, string> >& dest) {
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
}


void parseAlphabets(char* fileName) {
	ifstream in(fileName);
	getline(in, alphabets);
	for (int i=0;i<256;i++)
		mapAlphabets[i]=-1;
	for (int i=0;i<alphabets.length();i++)
		mapAlphabets[alphabets[i]] = i;
}


void parseScoreMatrix(char* fileName) {
	ifstream in(fileName);
	int tmp;
	for (int i=0;i<alphabets.length();i++) {
		for (int j=0;j<alphabets.length();j++) {
			cin >> tmp;
			score[i][j] = tmp;
		}
	}
}


AlignmentResult getAlignmentResult(vector < vector<int> >& dp, vector < vector<int> >& prev, int row, int col, string a, string b) {
	AlignmentResult res;
	res.sequence[0] = "";
	res.sequence[1] = "";
	res.score = dp[row][col];
	while(prev[row][col] != -1) {
		int r = prev[row][col]/(b.size()+1);
		int c = prev[row][col]%(b.size()+1);
		if (r == row-1 && c == col-1) {
			res.sequence[0].push_back(a[row]);
			res.sequence[1].push_back(b[col]);
		}
		else if (r == row-1) {
			res.sequence[0].push_back('.');
			res.sequence[1].push_back(b[col]);
		}
		else {
			res.sequence[0].push_back(a[row]);
			res.sequence[1].push_back('.');
		}
		row = r;
		col = c;
	}
	res.startPos[0] = row;
	res.startPos[1] = col;
	return res;
}


void intializeTables(vector < vector<int> >& dp, vector < vector<int> >& prev, int row, int col) {
	for (int i=0;i<row;i++) {
		dp.push_back(vector<int>(col));
		prev.push_back(vector<int>(col));
		dp[i][0] = dp[i-1][0] + gapPenality;
		prev[i][0] = -1;
	}
	
	for (int i=1;i<col;i++) {
		dp[0][i] = 0;
		prev[0][i] = -1;
	} 
}

void getMaxNeighbor(int &p, int& q, int& val, int i, int j, string& a, string& b, vector< vector<int> >& dp) {
	p=i-1;
	q=j-1;
	val = dp[p][q]+score[mapAlphabets[a[p]]][mapAlphabets[b[q]]];
	if (val < dp[i-1][j]+gapPenality) {
		val = dp[i-1][j]+gapPenality;
		p = i-1;
	}
	
	if (val < dp[i][j-1]+gapPenality) {
		val = dp[i][j-1]+gapPenality;
		q = j-1;
	}
}

AlignmentResult getGlobalAlignment(string a, string b) {
	int row = a.size()+1;
	int col = b.size()+1;
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
	getAlignmentResult(dp, prev, row-1, col-1, a, b);
}



AlignmentResult getLocalAlignment(string a, string b) {
	int row = a.size()+1;
	int col = b.size()+1;
	vector < vector<int> >dp;
	vector < vector<int> > prev;
	intializeTables(dp, prev, row, col);
	for (int i=1;i<row;i++) {
		for (int j=1;j<col;j++) {
			int p, q, val;
			getMaxNeighbor(p, q, val, i, j, a, b, dp);
			
			dp[i][j] = val;
			if (val <= 0) {
				prev[i][j] = -1;
			}
			else {
				prev[i][j] = p*col+q;
			}
		}
	}
	getAlignmentResult(dp, prev, row-1, col-1, a, b);
}

AlignmentResult getEndSpaceAlignment(string a, string b) {
	int row = a.size()+1;
	int col = b.size()+1;
	vector < vector<int> >dp;
	vector < vector<int> > prev;
	
	intializeTables(dp, prev, row, col);
	
	for (int i=1;i<row;i++) {
		for (int j=1;j<col;j++) {
			int p, q, val;
			getMaxNeighbor(p, q, val, i, j, a, b, dp);
			
			dp[i][j] = val;
			prev[i][j] = p*col+q;
		}
	}
	getAlignmentResult(dp, prev, row-1, col-1, a, b);
}

void topKalignments(int alignmentMethod) {
	map<int, vector<AlignmentResult> >mp;
	map<int, vector<int> >execMap;
	for (int i=0;i<queries.size();i++) {
		int start_s=clock();
		for (int j=0;j<database.size();j++) {
			AlignmentResult tmp;
			if (alignmentMethod == 1)
				tmp = getGlobalAlignment(queries[i].second, database[i].second);
			else if (alignmentMethod == 2)
				tmp = getLocalAlignment(queries[i].second, database[i].second);
			else
				tmp = getEndSpaceAlignment(queries[i].second, database[i].second);
			tmp.ids[1] = database[i].first;
			tmp.ids[0] = queries[i].first;
			mp[tmp.score].push_back(tmp);
		}
		start_s = (clock() - start_s)/double(CLOCKS_PER_SEC)*1000;
		execMap[queries[i].second.length()].push_back(start_s);
	}
	int cnt = 0;
	for (map<int, vector<AlignmentResult> >::iterator itr=mp.begin(); itr!= mp.end() && cnt < neighbors;itr++) {
		for (int i=0;i<itr->second.size() && cnt < neighbors; i++) {
			cnt++;
			AlignmentResult tmp = itr->second[i];
			cout<<"Score = "<<tmp.score<<endl;
			cout<<tmp.ids[0]<<" "<<tmp.startPos[0]<<" "<<tmp.sequence[0]<<endl;
			cout<<tmp.ids[1]<<" "<<tmp.startPos[1]<<" "<<tmp.sequence[1]<<endl;
		}
	}
}


int main(int argc, char *argv[]) {
	if (argc != 8) {
		cout<<"Invalid input"<<endl;
		return 0;
	}
	
	readSequencesFromFile(argv[2], queries);
	readSequencesFromFile(argv[3], database);
	parseAlphabets(argv[4]);
	parseScoreMatrix(argv[5]);
	neighbors = atoi(argv[6]);
	gapPenality = atoi(argv[7]);
	topKalignments(atoi(argv[0]));
    return 0;
}    
