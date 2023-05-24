#include <cstdio>
#include <algorithm>
#include <iostream>
using namespace std;
const int B=2;

int main() {
	int pf,qf,pt,qt,M;
	scanf("%d%d%d%d%d",&pf,&qf,&pt,&qt,&M);
	pair<int,int> A0=make_pair(pf,qf);
	pair<int,int> B0=make_pair(pt,qt);
	int ans = 0;
	
	for(int _pf=B;_pf>=-B;_pf--)
		for(int _qf=B;_qf>=-B;_qf--)
			for(int _pt=B;_pt>=-B;_pt--)
				for(int _qt=B;_qt>=-B;_qt--) {
					pair<int,int> A=make_pair(_pf,_qf);
					pair<int,int> B=make_pair(_pt,_qt);
					if(A<A0) continue;
					if(A==A0 && B<B0) continue;
					if(A>make_pair(1,0))continue;
					if(B<make_pair(-1,0))continue;
					if(A==A0 && B==B0) ans+=M; else if(A==make_pair(1,0)) ans+=5; else ans+=25;
					cout << A.first << "," << A.second << " " << B.first << "," << B.second << endl;
				}
	printf("%d\n",ans);
	return 0;
}