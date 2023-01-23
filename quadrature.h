struct Quadrature{
	int n;
	virtual double next() = 0;
};
template<class T>
struct Trapzd : Quadrature {
	double a,b,s;
	T &func;
	Trapzd() {};
	Trapzd(T &funcc, const double aa, const double bb) :
		func(funcc), a(aa), b(bb) {n=0;}
	double next() {
		double x,tnm,sum,del;
		int it,j;
		n++;
		if (n == 1) {
			return (s=0.5*(b-a)*(func(a)+func(b)));
		} else {
			for (it=1,j=1;j<n-1;j++) it <<= 1;
			tnm=it;
			del=(b-a)/tnm;
			x=a+0.5*del;
			for (sum=0.0,j=0;j<it;j++,x+=del) sum += func(x);
			s=0.5*(s+(b-a)*sum/tnm);
			return s;
		}
	}
};
template<class T>
double qtrap(T &func, const double a, const double b, const double eps=1.0e-7) {
	const int JMAX=20;
	double s,olds=0.0;
	Trapzd<T> t(func,a,b);
	for (int j=0;j<JMAX;j++) {
		s=t.next();
		if (j > 5)
			if (abs(s-olds) < eps*abs(olds) ||
				(s == 0.0 && olds == 0.0)) return s;
		olds=s;
	}
	throw("Too many steps in routine qtrap");
}
template<class T>
double qsimp(T &func, const double a, const double b, const double eps=1.0e-7) {
	const int JMAX=20;
	double s,st,ost=0.0,os=0.0;
	Trapzd<T> t(func,a,b);
	for (int j=0;j<JMAX;j++) {
		st=t.next();
		s=(4.0*st-ost)/3.0;
		if (j > 5)
			if (abs(s-os) < eps*abs(os) ||
				(s == 0.0 && os == 0.0)) return s;
		os=s;
		ost=st;
	}
	throw("Too many steps in routine qsimp");
}
template <class T>
struct Midpnt : Quadrature {
	double a,b,s;
	T &funk;
	Midpnt(T &funcc, const double aa, const double bb) :
		funk(funcc), a(aa), b(bb) {n=0;}
	double next(){
		int it,j;
		double x,tnm,sum,del,ddel;
		n++;
		if (n == 1) {
			return (s=(b-a)*func(0.5*(a+b)));
		} else {
			for(it=1,j=1;j<n-1;j++) it *= 3;
			tnm=it;
			del=(b-a)/(3.0*tnm);
			ddel=del+del;
			x=a+0.5*del;
			sum=0.0;
			for (j=0;j<it;j++) {
				sum += func(x);
				x += ddel;
				sum += func(x);
				x += del;
			}
			s=(s+(b-a)*sum/tnm)/3.0;
			return s;
		}
	}
	virtual double func(const double x) {return funk(x);}
};
template <class T>
struct Midinf : Midpnt<T>{
	double func(const double x) {
		return Midpnt<T>::funk(1.0/x)/(x*x);
	}
	Midinf(T &funcc, const double aa, const double bb) :
		Midpnt<T>(funcc, aa, bb) {
		Midpnt<T>::a=1.0/bb;
		Midpnt<T>::b=1.0/aa;
	}
};
template <class T>
struct Midsql : Midpnt<T>{
	double aorig;
	double func(const double x) {
		return 2.0*x*Midpnt<T>::funk(aorig+x*x);
	}
	Midsql(T &funcc, const double aa, const double bb) :
		Midpnt<T>(funcc, aa, bb), aorig(aa) {
		Midpnt<T>::a=0;
		Midpnt<T>::b=sqrt(bb-aa);
	}
};
template <class T>
struct Midsqu : Midpnt<T>{
	double borig;
	double func(const double x) {
		return 2.0*x*Midpnt<T>::funk(borig-x*x);
	}
	Midsqu(T &funcc, const double aa, const double bb) :
		Midpnt<T>(funcc, aa, bb), borig(bb) {
		Midpnt<T>::a=0;
		Midpnt<T>::b=sqrt(bb-aa);
	}
};
template <class T>
struct Midexp : Midpnt<T>{
	double  func(const double x) {
		return Midpnt<T>::funk(-log(x))/x;
	}
	Midexp(T &funcc, const double aa, const double  bb) :
		Midpnt<T>(funcc, aa, bb) {
		Midpnt<T>::a=0.0;
		Midpnt<T>::b=exp(-aa);
	}
};