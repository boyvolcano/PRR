#include <map>
#include <cassert>
#include <cmath>
#include <vector>
#include <string>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <algorithm>
using namespace std;
int B = 8;

struct vec6 {

	int deg_n,deg_ln;
	int deg_a,deg_la;
	int deg_v,deg_lv;

	vec6() = default;

	vec6(const vec6 &a) = default;

	vec6 operator + (const vec6 & a) const {
		return {deg_n + a.deg_n, deg_ln + a.deg_ln, deg_a + a.deg_a, deg_la + a.deg_la, deg_v + a.deg_v, deg_lv + a.deg_lv};
	}

	vec6 operator - (const vec6 & a) const {
		return {deg_n - a.deg_n, deg_ln - a.deg_ln, deg_a - a.deg_a, deg_la - a.deg_la, deg_v - a.deg_v, deg_lv - a.deg_lv};		
	}

	int operator < (const vec6 & a) const {
		if(deg_n != a.deg_n) return deg_n < a.deg_n;
		if(deg_ln != a.deg_ln) return deg_ln < a.deg_ln;
		if(deg_v != a.deg_v) return deg_v < a.deg_v;
		if(deg_lv != a.deg_lv) return deg_lv < a.deg_lv;
		if(deg_a != a.deg_a) return deg_a < a.deg_a;
		if(deg_la != a.deg_la) return deg_la < a.deg_la;
		return 0;
	}

	vec6 substitute() const {
		return {0, 0, deg_a, deg_la, deg_n, deg_ln};
	}

	void print() const {
		if(deg_n) printf("n^(%d)*",deg_n);
		if(deg_ln) printf("(ln n)^(%d)*",deg_ln);
		if(deg_v) printf("v^(%d)*",deg_v);
		if(deg_lv) printf("(ln v)^(%d)*",deg_lv);
		if(deg_a) printf("a^(%d)*",deg_a);
		if(deg_la) printf("(ln a)^(%d)*",deg_la);
	}

	int get_dependency() const {
		int x = 0;
		if(deg_n) x |= 1;
		if(deg_ln) x |= 2;
		if(deg_a) x |= 4;
		if(deg_la) x |= 8;
		if(deg_v) x |= 16;
		if(deg_lv) x |= 32;
		return x;
	}

	int independent(int msk) const {
		return (get_dependency() & msk) == 0;
	} 
};

struct poly{
	map<vec6, long double> coef;

	poly() {
		coef.clear();
	}

	~poly() {
		coef.clear();
	}

	poly (const poly & a) {
		for(const auto & u: a.coef) {
			coef[u.first] = u.second;
		}
	}

	poly operator +(const poly &a) const {
		poly res;
		for(const auto & u: a.coef)
			res.coef[u.first] += u.second;
		for(const auto & u: coef)
			res.coef[u.first] += u.second;
		return res;
	}

	poly operator -(const poly &a) const {
		poly res;
		for(const auto & u: a.coef)
			res.coef[u.first] -= u.second;
		for(const auto & u: coef)
			res.coef[u.first] += u.second;
		return res;
	}

	poly operator /(const poly &a) const {
		poly res;
		assert(a.check_monomial());
		for(const auto & u: a.coef) {
			for(const auto & v: coef)
				res.coef[v.first - u.first] += v.second / u.second;
			break;
		}
		return res;
	}

	poly operator *(const poly &a) const {
		poly res;
		for(const auto & u: a.coef) {
			for(const auto & v: coef)
				res.coef[v.first + u.first] += v.second * u.second;
		}
		return res;
	}

	poly operator *(const long double &a) const {
		poly res;
		for(const auto & v: coef)
			res.coef[v.first] += v.second * a;
		return res;

	}

	void add(vec6 x, long double y) {
		coef[x] += y;
	}

	poly substitute() const {
		poly res;
		for(const auto& u: coef) 
			res.coef[u.first.substitute()] += u.second;
		return res;
	}

	poly scale(long double lambda) const {
		poly res;
		long double llambda = log(lambda);
		for(const auto& u: coef) { 
			long double c = u.second;
			vec6 p = u.first;
			if(u.first.deg_n > 0)
				for(int j = 0; j < u.first.deg_n; j ++) c *= lambda;
			if(u.first.deg_n < 0)
				for(int j = 0; j < -u.first.deg_n; j ++) c /= lambda;
			for(int k = 0; k <= u.first.deg_ln; k ++) {
				if(k) c = c / k * (u.first.deg_ln - k + 1);
				long double d = c;
				for(int l = 0; l < k; l ++) d = d * llambda;
				res.coef[{p.deg_n, p.deg_ln - k, p.deg_a, p.deg_la, p.deg_v, p.deg_lv}] += d;
			}
		}
		return res;
	}

	void print() const {
		for(const auto &u: coef) {
			if(u.second) {
				u.first.print();
				printf(" = %Lf\n", u.second);
			}
		}
	}

	int independent(int msk) const {
		for(const auto& u: coef) 
			if(!u.first.independent(msk))
				return 0;
		return 1;
	}

	bool check_monomial() const {
		return coef.size() == 1;
	}

	bool is_empty() const {
		return coef.size() == 0;
	}

	pair<vec6, long double> sig_term() const {
		return *(--coef.end());
	}

	void simplified() {
		map<vec6, long double> tmp = coef;
		coef.clear();
		for(const auto& x: tmp) {
			if(fabsl(x.second) >= 1e-6) {
				coef[x.first] = x.second;
			}
		}
	}
};

struct canonical {
	vector< pair<long double, poly> > expr;

	canonical() {
		expr.clear();
	}

	canonical(const canonical& a): expr(a.expr) {}

	~canonical() { expr.clear(); }

	canonical operator + (const canonical& a) const {
		canonical res;
		for(const auto& v: expr)
			res.expr.emplace_back(v);
		for(const auto& v: a.expr)
			res.expr.emplace_back(v);
		return res;
	}
	canonical operator *(const double& w) const {
		canonical res;
		for(const auto& u: expr) 
			res.expr.emplace_back(u.first * w, u.second);
		return res;
	}

	void print() const {
		for(const auto& term: expr) {
			if(term.first > 0.0000001) {
				printf("%Lf\n", term.first);
				term.second.print();
				puts("");
			}
		}
	}
};


struct distribution {
	string name;

	distribution() {}

	distribution(const distribution& a): name(a.name) {} 

	distribution(const string& a): name(a) {}

	~distribution() { name.clear(); }

	virtual distribution* getcopy() const {return nullptr;}
};

struct discrete: public distribution {
	vector< pair<long double, poly> > expr;

	discrete(): distribution("discrete") {
		expr.clear();
	}

	discrete(const discrete& a): expr(a.expr), distribution("discrete") {}

	~discrete() {
		expr.clear();
	}

	virtual distribution* getcopy() const {
		return new discrete(*this);
	}
};


struct uniform: public distribution {
	uniform(): distribution("uniform") {}
	uniform(const uniform& a): distribution("uniform") {}


	virtual distribution* getcopy() const {
		return new uniform;
	}
};


struct muniform: public distribution {
	muniform(): distribution("muniform") {}
	muniform(const muniform& a): distribution("muniform") {}

	virtual distribution* getcopy() const {
		return new muniform;
	}

};


struct rec_body {
	poly pretime;
	distribution* v_dist;
	int rec_type, b, c;

	rec_body() {
		pretime = poly();
		rec_type = b = c = 0;
		v_dist = nullptr;
	}

	rec_body(const rec_body& a):pretime(a.pretime) {
		rec_type = a.rec_type;
		b = a.b;
		c = a.c;
		v_dist = a.v_dist->getcopy();
	}

	~rec_body() {
		delete v_dist;
	}

	void read() {
		pretime = poly();
		rec_type = b = c = 0;
	}

};

poly monomial(vec6 p, long double u) {
	poly res;
	res.add(p,u);
	return res;
}


poly replace_const(const poly &a, long double u) {
	
	poly res;
	static long double cur_n[10], cur_ln[10];
	long double lu = logl(u);
	cur_n[0] = cur_ln[0] = 1;
	for(int i = 1; i < 10; i ++)
		cur_n[i] = cur_n[i-1] * u;
	for(int i = 1; i < 10; i ++)
		cur_ln[i] = cur_ln[i-1] * lu;

	for(const auto& term: a.coef) {
		vec6 p(term.first);
		p.deg_n = p.deg_ln = 0;
		res = res + monomial(p,term.second * cur_n[term.first.deg_n] * cur_ln[term.first.deg_ln]);
	}
	
	return res;
}

poly replace(const poly &a, const poly& u) {
	
	assert(u.independent(62));
	assert(a.independent(2));

	poly res;
	static poly cur[10];
	cur[0] = monomial({0,0,0,0,0,0},1);
	for(int i = 1; i < 5; i ++)
		cur[i] = cur[i-1] * u;

	for(const auto& term: a.coef) {
		vec6 p(term.first);
		p.deg_n = p.deg_ln = 0;
		res = res + monomial(p,term.second) * cur[term.first.deg_n];
	}
	
	return res;
}

canonical discrete_approx(const rec_body &rec, const discrete& prb, const poly& t, const poly& f) {
	//puts("here");
	canonical res;

	for(const auto& u: prb.expr) {

		res.expr.emplace_back(u.first, t * (rec.pretime + replace(f, u.second) - f));
	}

	return res;
}

bool doable_single(const poly& w) {
	if(!w.check_monomial()) return 0;
	vec6 mono_term = w.sig_term().first;
	if(mono_term.deg_v)
		return 1; 
	if(mono_term.deg_lv)
		return 1; 
	return 0;
}

canonical uniform_approx(const rec_body& rec, const uniform& u, const poly& t, const poly& f) {
	canonical res;
	if(rec.rec_type == 1 || rec.rec_type == 2) {
		poly w = t * f.substitute();
		w.simplified();
		if(doable_single(w)) {
			pair<vec6, long double> vterm = w.sig_term();
			vec6 mono_term = vterm.first;
			if(mono_term.deg_v) {
				mono_term.deg_n += mono_term.deg_v - 1; mono_term.deg_v = 0;
				mono_term.deg_ln += mono_term.deg_lv; mono_term.deg_lv = 0;
				vec6 nwterm = mono_term;
				nwterm.deg_n ++;
				//printf("nterm = ");nwterm.print();puts("");
				//printf("mterm = ");mono_term.print();puts("");
				res.expr.emplace_back(1 / vterm.second, 
					t * rec.pretime - t * f +
					monomial(nwterm, vterm.second) + 
					monomial({0, 1, 0, 0, 0, 0}, -mono_term.deg_n-1 + max(0, -mono_term.deg_ln)) +
					monomial({0, 0, 0, 1, 0, 0}, -mono_term.deg_a + max(0, -mono_term.deg_la)));
			} else if(mono_term.deg_lv){
				mono_term.deg_n += mono_term.deg_v; mono_term.deg_v = 0;
				mono_term.deg_ln += mono_term.deg_lv - 1; mono_term.deg_lv = 0;
				vec6 nwterm = mono_term;
				nwterm.deg_ln ++;
				//printf("nterm = ");nwterm.print();puts("");
				//printf("mterm = ");mono_term.print();puts("");
				res.expr.emplace_back(1 / vterm.second, 
					t * rec.pretime - t * f +
					monomial(nwterm, vterm.second) + 
					monomial({0, 1, 0, 0, 0, 0}, -mono_term.deg_n + max(0, -mono_term.deg_ln)) +
					monomial({0, 0, 0, 1, 0, 0}, -mono_term.deg_a + max(0, -mono_term.deg_la)));
			}
		} else {
			assert(0);
		}
	} else {
		poly w = t * f.substitute();
		w.simplified();
		static int B = 4;
		if(w.check_monomial()) {
			pair<vec6, long double> vterm = w.sig_term();
			vec6 mono_term = vterm.first;
			if(mono_term.deg_v == 1 && mono_term.deg_lv == 0){
				res.expr.emplace_back(1, t * rec.pretime - t * vterm.second);
			} else {
				//f.scale(0.5).print(); puts("");
				for(int i = 0; i < B; i ++)
					res.expr.emplace_back(1./(2*B),t*(rec.pretime+f.scale((i+1.)/(2*B))+f.scale((2*B-i-1.)/(2*B))-f));
				for(int i = B; i < 2 * B; i ++)
					res.expr.emplace_back(1./(2*B),t*(rec.pretime+f.scale(i/(2.*B))+f.scale((2*B-i)/(2.*B))-f));
			}
		}
	}
	//res.print();
	return res;
}

canonical muniform_approx(const rec_body& rec, const muniform& m, const poly& t, const poly& f) {
	canonical res;
	if(rec.rec_type == 1 || rec.rec_type == 2) {
		poly w = t * f.substitute();
		w.simplified();
		if(doable_single(w)) {
			pair<vec6, long double> vterm = w.sig_term();
			vec6 mono_term = vterm.first;
			//w.print();
			if(mono_term.deg_v) {
				mono_term.deg_n += mono_term.deg_v - 1; mono_term.deg_v = 0;
				mono_term.deg_ln += mono_term.deg_lv; mono_term.deg_lv = 0;
				vec6 nwterm = mono_term;
				nwterm.deg_n ++;
				//printf("nterm = ");nwterm.print();puts("");
				//printf("mterm = ");mono_term.print();puts("");
				res.expr.emplace_back(2 / vterm.second, 
					t * rec.pretime - t * f +
					monomial(nwterm, vterm.second) + 
					monomial({0, 1, 0, 0, 0, 0}, -mono_term.deg_n-1 + max(0, -mono_term.deg_ln)) +
					monomial({0, 0, 0, 1, 0, 0}, -mono_term.deg_a + max(0, -mono_term.deg_la)));
			} else if(mono_term.deg_lv){
				mono_term.deg_n += mono_term.deg_v; mono_term.deg_v = 0;
				mono_term.deg_ln += mono_term.deg_lv - 1; mono_term.deg_lv = 0;
				vec6 nwterm = mono_term;
				nwterm.deg_ln ++;
//				printf("nterm = ");nwterm.print();puts("");
//				printf("mterm = ");mono_term.print();puts("");
				res.expr.emplace_back(2 / vterm.second, 
					t * rec.pretime - t * f +
					monomial(nwterm, vterm.second) + 
					monomial({0, 1, 0, 0, 0, 0}, -mono_term.deg_n + max(0, -mono_term.deg_ln)) +
					monomial({0, 0, 0, 1, 0, 0}, -mono_term.deg_a + max(0, -mono_term.deg_la)));
			}
		} else {
			assert(0);
		}
	} else {
		assert(0);
	}
	//res.print();
	return res;
}

rec_body read_body() {
	rec_body res;

	int a,b,c,d,e,f;long double g; string name;

	cin>>a>>b>>c>>d>>e>>f>>g;
	res.pretime = monomial({a,b,c,d,e,f},g);

	cin>>name;
	if(name == "M") res.v_dist = new muniform;
	if(name == "D"){
		discrete* p = new discrete;
		res.v_dist = p;
		int m;
		cin >> m;
		for(int i = 0; i < m; i ++) {
			long double prb;
			cin >> prb >> b;
			poly q = monomial({1,0,0,0,0,0},1) + monomial({0,0,0,0,0,0},b);
			p->expr.emplace_back(prb, q);
		}
	}
	if(name == "U") res.v_dist = new uniform;
	cin >> res.rec_type >> res.b >> res.c;

	return res;
}


struct file {

	int n_branch;
	vector< pair<long double, rec_body> > prog;
	poly tm_t, tm_f;
	int uf,vf,ut,vt;

	file() = default;

	void read() {
		cin >> n_branch;
		prog.clear();
		for(int i = 0; i < n_branch; i ++) {
			long double prb;
			cin >> prb;
			//cout << prb << endl;
			rec_body body = read_body();
			prog.emplace_back(prb, body);
		}
		int a,b,c,d,e,f;
		cin>>a>>b>>c>>d>>e>>f; uf=a;vf=b;//tm_f = monomial({0,0,c,d,0,0},1);
		cin>>a>>b>>c>>d>>e>>f; ut=a;vt=b;//tm_t = monomial({0,0,c,d,0,0},1);
	}
};

canonical file_approx(const file& fi, const poly& t, const poly& f) {
	//t.print();
	//f.print();
	canonical res;
	//puts("here2");
	for(const auto& term: fi.prog) {
		if(term.second.v_dist->name == "discrete") {
			canonical w = discrete_approx(term.second, *static_cast<discrete*>(term.second.v_dist), t, f);
			res = res + w * term.first;
		}
		if(term.second.v_dist->name == "uniform") {
			canonical w = uniform_approx(term.second, *static_cast<uniform*>(term.second.v_dist), t, f);
			res = res + w * term.first;
		}
		if(term.second.v_dist->name == "muniform") {
			canonical w = muniform_approx(term.second, *static_cast<muniform*>(term.second.v_dist), t, f);
			res = res + w * term.first;
		}
	}
	return res;
}

int check_ind_n(const canonical& res) {
	long double mu = 0;
	for(auto p: res.expr) {
		assert(p.second.independent(51));
		if(p.second.is_empty()) {
			mu += p.first;
			continue;
		}
		const auto& tm = p.second.sig_term();
		if(tm.first.independent(63)) {
			mu += p.first * expl(tm.second);
			continue;
		}
		if(p.first < 0)
			continue;
		if(tm.second > 0)
			return 0;
	}
//	printf("mu = %Lf\n", mu);
	return mu <= 1;
}

int check(const canonical& res) {
	for(int n = 2; n <= 10; n ++) {
		canonical res2 = res;
		for(auto &q: res2.expr){
			q.second = replace_const(q.second, n);
			q.second.simplified();
		}
		if(!check_ind_n(res2))
			return 0;
	}
	return 1;
}


int main () {
	//freopen("manual.txt", "r", stdin);
	file f;
	f.read();
	//puts("done");

	for(int _pt=2;_pt>=-2;_pt--)
		for(int _qt=2;_qt>=-2;_qt--)
			for(int _pf=-2;_pf<=2;_pf++)
				for(int _qf=-2;_qf<=2;_qf++) {
					pair<int,int> A=make_pair(_pf,_qf);
					pair<int,int> B=make_pair(_pt,_qt);
					if(A>make_pair(1,0))continue;
					if(B<make_pair(-1,0))continue;
					f.tm_f = monomial({f.uf,f.vf,_pf,_qf,0,0}, 1);
					f.tm_t = monomial({f.ut,f.vt,_pt,_qt,0,0}, 1);
					int L = 16;
					//f.tm_f.print();
					//f.tm_t.print();
					if(A==make_pair(1,0))L=1; else L=16;
					for(int j = 1; j <= 16; j *= 2)
					for(int i = 1; i <= L; i *= 2)
						 {
							canonical res = file_approx(f, f.tm_t * (1.0 / j), f.tm_f * (i * 0.5));
							//res.print();
							if(check(res)) {
								printf("%d %d\n",i,j);
								((f.tm_t * (1.0 / j)) * (f.tm_f * (i * 0.5))).print();
								((f.tm_t * (1.0 / j)) * monomial({f.uf, f.vf, 1,0,0,0},1)).print();
								return 0;
							}
						}
				}
}
