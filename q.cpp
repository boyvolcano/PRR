#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <random>
using namespace std;

default_random_engine e(time(0));
uniform_int_distribution<int> u(1,15);

int main() {
	int n; scanf("%d",&n);
	long double t=0;
	while(n--){
		t+=u(e)*0.001;
	}
	cout<<t<<endl;
	return 0;
}
