#include <cstdlib>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include <math.h>
#include <vector>
#include <string>
#include <random>
#include <algorithm>

using namespace std;

//generates a pseudo-random double between 0.0 and 0.999...
double randdouble()
{
    return rand()/(double(RAND_MAX)+1);
}

// generates realization of binomial point process
vector<vector<double>> binomial_process(int n) {
	vector<vector<double>> pattern(n,vector<double>(2));
	for (int i=0;i<n;++i) {
		pattern[i][0]=randdouble();
		pattern[i][1]=randdouble();
	}
	return pattern;
}

// check for overlap between circles of radius r
bool check_overlap(vector<double> xk_trial, vector<vector<double>> pattern, int k, int n, double r) {
	bool overlap=false;
	
	for (int i=0;i<n;++i) {
		if (i!=k) {
			if (pow(pattern[i][0]-xk_trial[0],2)+pow(pattern[i][1]-xk_trial[1],2)<pow(2.*r,2)) {
				overlap=true;
				break;
			}
		}
	}
	return overlap;
}

// write output file of marked point pattern
int write_output(int n,vector<vector<double>> pattern,vector<vector<double>> marks,string file_name) {
	ofstream myfile;
	myfile.open (file_name);
	for (int i=0;i<n;++i) {
		myfile << pattern[i][0] << "," << pattern[i][1] << "," << marks[i][0] << "," << marks[i][1] << "," << marks[i][2] << endl;
	}
	myfile.close();
	return 0;
}

// write output file of sRx
int write_output2(int n,vector<int> sRx,string file_name) {
	ofstream myfile;
	myfile.open (file_name);
	for (int i=0;i<n;++i) {
		myfile << sRx[i] << endl;
	}
	myfile.close();
	return 0;
}


///////////////////////////////////////// utility ////////////////////////////////////////////
//
// split a string
vector<string> split(const string &s, char delim) {
    stringstream ss(s);
    string item;
    vector<string> tokens;
    while (getline(ss, item, delim)) {
        tokens.push_back(item);
    }
    return tokens;
}
//
// sgn function
int sign(double val) {
	int sgn=0;
	if (val>0) {
		return 1;
	}
	else if (val<0) {
		return -1;
	}
	return sgn;
}
//
// linspace function
vector<double> linspace(double a, double b, int n) {
    vector<double> array;
    double dt=(b-a)/(n-1);	
	for (int i=0;i<n;++i) {
		array.push_back(a+i*dt);
	}
    return array;
}

const double pi=4*atan(1.);
//
// refine indices to account for the fact that a -1 index does not point
// to the last element of a vector -- contrarily to Python implementation
int re_j(int j,int size) {
	if (j==0) {
		return j+size;
	}
	else {
		return j;
	}
}
//
// factorial function
int factorial(int n) {
  return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}
//
// binomial coefficient
int binom(int n, int k) {
    int res = 1;
    // Since C(n, k) = C(n, n-k)
    if ( k > n - k )
        k = n - k;
    // Calculate value of [n * (n-1) *---* (n-k+1)] / [k * (k-1) *----* 1]
    for (int i = 0; i < k; ++i) {
        res *= (n - i);
        res /= (i + 1);
    }
    return res;
}
//
// get non-repeated permutations
vector<vector<int>> get_perms(vector<int> inds) {
	vector<vector<int>> perms;
	sort(inds.begin(),inds.end());
	do {
		perms.push_back(inds);
	} while (next_permutation(inds.begin(),inds.end()));
	return perms;
}
//
// compute the "fully" symmetric tensor product a^{\otimes l}\otimes b^{\otimes m}\otimes c^{\otimes n}
double get_abc_sym(vector<double> a,vector<double> b,vector<double> c,int l, int m, int n, vector<int> inds) {
	double comp=0.;
	int order=inds.size();
	int n1=0,n2=0;
	for (int i=0;i<order;++i) {
		if (inds[i]==1) {++n1;}
		else if (inds[i]==2) {++n2;}
	}
	double p=double(factorial(order))/factorial(n1)/factorial(n2);
	vector<vector<int>> perms=get_perms(inds);
	for (size_t i=0;i<perms.size();++i) {
		vector<int> ind=perms[i];
		double comp_ind=1.;
		for (int k_ind=0;k_ind<l;++k_ind) {
			comp_ind*=a[ind[k_ind]-1];
		}
		for (int k_ind=l;k_ind<l+m;++k_ind) {
			comp_ind*=b[ind[k_ind]-1];
		}
		for (int k_ind=l+m;k_ind<l+m+n;++k_ind) {
			comp_ind*=c[ind[k_ind]-1];
		}
		comp+=comp_ind;
	}
	comp/=p;
	return comp;
}
//
// compute integral i_st_k for w2_rs
double int_k(int s,int t,double dt) {
	double ival;
	if (t==0) {
		if (s==0) {
			ival=dt;
		}
		else if (s==1) {
			ival=sin(dt);
		}
		else if (s==2) {
			ival=.5*(dt+sin(dt)*cos(dt));
		}
		else if (s==3) {
			ival=1./12.*(9.*sin(dt)+sin(3.*dt));
		}
		else if (s==4) {
			ival=1./32.*(12.*dt+8.*sin(2.*dt)+sin(4.*dt));
		}
		else if (s==5) {
			ival=1./8.*(5.*sin(dt)+5./6.*sin(3.*dt)+1./10.*sin(5.*dt));
		}
		else if (s==6) {
			ival=1./192.*(60.*dt+45.*sin(2.*dt)+9.*sin(4.*dt)+sin(6.*dt));
		}
	}
	else if (t==1) {
		if (s>=1) {
			ival=1./s*(1.-pow(cos(dt),s));
		}
	}
	else if (t==2) {
		if (s==2) {
			ival=1./2.*(dt-sin(dt)*cos(dt));
		}
		else if (s==3) {
			ival=1/3.*pow(sin(dt),3);
		}
		else if (s==4) {
			ival=1./32.*(4.*dt-sin(4.*dt));
		}
	}
	else if (t==3) {
		if (s==3) {
			ival=1./12*(cos(3.*dt)-9.*sin(dt)+8.);
		}
	}
	return ival;
}

///////////////////////////////////////// bnd_resovle ////////////////////////////////////////////


//////////////////////////////// minkowski characterization
//
// data structure for Minkowski tensor of order r+s and type nu
struct W_rs_nu {
	vector<double> comps;
	vector<vector<int>> inds;
};
//
// compute Minkowski tensor W_r_0 of a polytope
W_rs_nu get_W0_rs_type_II_k(int r, int s, vector<double> vk, vector<double> vk1, vector<double> nk, double Lk) {
	vector<double> W0_rs;
	vector<vector<int>> W0_rs_inds;
	for (int dd=0;dd<(r+s)/2+1;++dd) {		
		for (int d=0;d<2;++d) {
			if (int(W0_rs.size())<=r+s) {
				vector<int> inds(r+s);
				for (int k=0;k<r+s-dd;++k) {inds[k]=d+1;}
				for (int k=r+s-dd;k<r+s;++k) {inds[k]=2-d;}
				W0_rs.push_back(0.);
				W0_rs_inds.push_back(inds);
				for (int i=0;i<r+2;++i) {
					for (int j=0;j<i+1;++j) {
						vector<int> inds_p1=inds;
						inds_p1.push_back(1);
						W0_rs[W0_rs.size()-1]+=1./(2.+r)*binom(r+1,i)*binom(i,j)*pow(-1,i-j)*Lk/(1.+i)*get_abc_sym(vk,vk1,vk,r-j+1,j,0,inds_p1)*nk[0];
						inds_p1.back()=2;
						W0_rs[W0_rs.size()-1]+=1./(2.+r)*binom(r+1,i)*binom(i,j)*pow(-1,i-j)*Lk/(1.+i)*get_abc_sym(vk,vk1,vk,r-j+1,j,0,inds_p1)*nk[1];
					}
				}
			}
		}
	}
	return {W0_rs,W0_rs_inds};
}
//
// compute Minkowski tensor W_rs_1 of a polytope
W_rs_nu get_W1_rs_type_II_k(int r, int s, vector<double> vk, vector<double> vk1, vector<double> nk, double Lk) {
	vector<double> W1_rs;
	vector<vector<int>> W1_rs_inds;
	for (int dd=0;dd<(r+s)/2+1;++dd) {		
		for (int d=0;d<2;++d) {
			if (int(W1_rs.size())<=r+s) {
				vector<int> inds(r+s);
				for (int k=0;k<r+s-dd;++k) {inds[k]=d+1;}
				for (int k=r+s-dd;k<r+s;++k) {inds[k]=2-d;}
				W1_rs.push_back(0.);
				W1_rs_inds.push_back(inds);
				for (int i=0;i<r+1;++i) {
					for (int j=0;j<i+1;++j) {
						W1_rs[W1_rs.size()-1]+=.5*(Lk/(1.+i))*binom(r,i)*binom(i,j)*pow(-1,i-j)*get_abc_sym(vk,vk1,nk,r-j,j,s,inds);						
					}
				}
			}
		}
	}
	return {W1_rs,W1_rs_inds};
}
//
// compute Minkowski tensor W_rs_2 of a polytope
W_rs_nu get_W2_rs_type_II_k(int r, int s, vector<double> vk, vector<double> vk1, vector<double> nk, double Lk, double dtk) {
	vector<double> W2_rs;
	vector<vector<int>> W2_rs_inds;
	for (int dd=0;dd<(r+s)/2+1;++dd) {		
		for (int d=0;d<2;++d) {
			if (int(W2_rs.size())<=r+s) {
				vector<int> inds(r+s);
				for (int k=0;k<r+s-dd;++k) {inds[k]=d+1;}
				for (int k=r+s-dd;k<r+s;++k) {inds[k]=2-d;}
				W2_rs.push_back(0.);
				W2_rs_inds.push_back(inds);
				for (int i=0;i<s+1;++i) {
					for (int j=0;j<i+1;++j) {
						W2_rs[W2_rs.size()-1]+=.5*pow(Lk,-i)*binom(s,i)*binom(i,j)*pow(-1,i-j)*get_abc_sym(vk,vk1,nk,i-j,r+j,s-i,inds)*int_k(s,i,dtk);
					}
				}
			}
		}
	}
	return {W2_rs,W2_rs_inds};
}

//////////////////////////////// read import_bnd_solve_data
// define data structure for data read through import_bnd_solve_data
struct bnd_solve_data {
	int ngrains;
	vector<int> igrain_ids;
	vector<vector<int>> lorder;
	vector<vector<double>> c_ij;
	vector<vector<vector<vector<vector<double>>>>> thsr_ij;
	vector<vector<double>> cmn_thsr;
	vector<vector<vector<int>>> cmn_types;
	vector<vector<vector<double>>> uij1,uij2,zij,xij;
};
// function to read in data from bnd_solve
bnd_solve_data import_bnd_solve_data(string fname_in) {
	ifstream file_in;
	file_in.open(fname_in);
	string kline;
	getline(file_in,kline);
	int ncells=atoi(kline.c_str());
	int nneighbs=0;
	double tiny=pow(10.,-12);
	vector<int> igrain_ids;
	vector<vector<int>> list_lorder;
	vector<vector<double>> list_c_ij;
	vector<vector<vector<vector<vector<double>>>>> list_thsr_ij;
	vector<vector<double>> list_cmn_thsr;
	vector<vector<vector<int>>> list_cmn_types;
	vector<vector<vector<double>>> list_uij1;
	vector<vector<vector<double>>> list_uij2;
	vector<vector<vector<double>>> list_zij;
	vector<vector<vector<double>>> list_xij;
	for (int k=1;k<ncells+1;++k) {
		getline(file_in,kline); // read k-th line
		vector<string> x=split(kline,','); // split line between ','
		igrain_ids.push_back(atoi(x[0].c_str())); // populate igrain_id
		nneighbs=atoi(x[1].c_str()); // read n_neighbs
		//
		vector<int> lorder;
		for (int kk=2;kk<nneighbs+3;++kk) {
			lorder.push_back(atoi(x[kk].c_str()));
		}
		list_lorder.push_back(lorder); // populate list_lorder
		//
		vector<double> c_ij;
		for (int kk=nneighbs+3;kk<2*nneighbs+3;++kk) {
			c_ij.push_back(atof(x[kk].c_str()));
		}
		list_c_ij.push_back(c_ij); // populate list_c_ij
		//
		int i=2*nneighbs+3;
		//
		vector<vector<vector<vector<double>>>> thsr_ij;
		for (int kk=0;kk<nneighbs;++kk) {
			int n_kk=atoi(x[i].c_str());
			vector<vector<vector<double>>> thsr_ij_kk;
			for (int jj=0;jj<n_kk;++jj) {
				i+=1;
				int n_jj=atoi(x[i].c_str());
				vector<vector<double>> thsr_ij_kk_jj(n_jj,vector<double>(2));
				for (int ii=0;ii<n_jj;++ii) {
					thsr_ij_kk_jj[ii][0]=atof(x[i+1].c_str())+tiny;
					thsr_ij_kk_jj[ii][1]=atof(x[i+2].c_str())-tiny;
					i+=2;
				}
				thsr_ij_kk.push_back(thsr_ij_kk_jj);
			}
			i+=1;
			thsr_ij.push_back(thsr_ij_kk);
		}
		list_thsr_ij.push_back(thsr_ij); // populate list_thsr_ij
		//		
		vector<double> cmn_thsr;
		for (int kk=i;kk<i+nneighbs;++kk) {
			cmn_thsr.push_back(atof(x[kk].c_str()));
		}
		list_cmn_thsr.push_back(cmn_thsr); // populate list_cmn_thsr		
		i+=nneighbs;	
		//
		vector<vector<int>> cmn_types(nneighbs,vector<int>(2));
		for (int kk=0;kk<nneighbs;++kk) {
			cmn_types[kk][0]=atoi(x[i+2*kk].c_str());
			cmn_types[kk][1]=atoi(x[i+2*kk+1].c_str());
		}
		list_cmn_types.push_back(cmn_types); // populate list_cmn_types
		i+=2*nneighbs;		
		//		
		vector<vector<double>> uij1(nneighbs,vector<double>(2));
		for (int kk=0;kk<nneighbs;++kk) {
			uij1[kk][0]=atof(x[i+2*kk].c_str());
			uij1[kk][1]=atof(x[i+2*kk+1].c_str());
		}
		list_uij1.push_back(uij1); // populate list_uij1
		i+=2*nneighbs;
		//
		vector<vector<double>> uij2(nneighbs,vector<double>(2));
		for (int kk=0;kk<nneighbs;++kk) {
			uij2[kk][0]=atof(x[i+2*kk].c_str());
			uij2[kk][1]=atof(x[i+2*kk+1].c_str());
		}
		list_uij2.push_back(uij2); // populate list_uij2
		i+=2*nneighbs;
		//
		vector<vector<double>> zij(nneighbs,vector<double>(3));
		for (int kk=0;kk<nneighbs;++kk) {
			zij[kk][0]=atof(x[i+3*kk].c_str());
			zij[kk][1]=atof(x[i+3*kk+1].c_str());
			zij[kk][2]=atof(x[i+3*kk+2].c_str());
		}
		list_zij.push_back(zij); // populate list_zij
		i+=3*nneighbs;
		//		
		vector<vector<double>> xij(nneighbs,vector<double>(2));
		for (int kk=0;kk<nneighbs;++kk) {	
			xij[kk][0]=atof(x[i+2*kk].c_str());
			xij[kk][1]=atof(x[i+2*kk+1].c_str());		
		}
		list_xij.push_back(xij); // populate list_xij
	}
	//
	file_in.close();
	//
	return {ncells,igrain_ids,list_lorder,list_c_ij,list_thsr_ij,list_cmn_thsr,list_cmn_types,list_uij1,list_uij2,list_zij,list_xij};
}	


//////////////////////////////// read marked point pattern
// define data structure for a marked point pattern (mpp)
struct marked_point_pattern {
	int ngrains;
	double lx,ly;
	vector<double> a,b,th;
	vector<vector<double>> x;
};
// function to read a marked point pattern file
marked_point_pattern read_mpp(string fname_in) {
	vector<double> a,b,th;
	vector<vector<double>> x;
	//
	ifstream file_in;
	file_in.open(fname_in);
	string kline;
	getline(file_in,kline);
	int ncells=atoi(kline.c_str());
	getline(file_in,kline);
	vector<string> x_str=split(kline,','); // split line between ','
	double lx=-atof(x_str[0].c_str());
	double ly=-atof(x_str[1].c_str());	
	getline(file_in,kline);
	x_str=split(kline,','); // split line between ','
	lx+=atof(x_str[0].c_str());
	ly+=atof(x_str[1].c_str());	
	//
	for (int k=0;k<ncells;++k) {
		getline(file_in,kline); // read k-th line
		x_str=split(kline,','); // split line between ','		
		vector<double> xk(2);
		xk[0]=atof(x_str[0].c_str());
		xk[1]=atof(x_str[1].c_str());
		x.push_back(xk);
		a.push_back(atof(x_str[2].c_str()));
		b.push_back(atof(x_str[2].c_str())*atof(x_str[3].c_str()));
		th.push_back(atof(x_str[4].c_str()));
	}
	//
	file_in.close();
	//
	return {ncells,lx,ly,a,b,th,x};
}


//////////////////////////////// read lists of neighbors from gtess
// function to read a list of neighbors obtained from gtess
vector<vector<int>> read_neighbs(string fname_in) {
	vector<vector<int>> neighbs;
	//
	ifstream file_in;
	file_in.open(fname_in);
	string kline;
	getline(file_in,kline);
	int ncells=atoi(kline.c_str());
	for (int k=0;k<ncells;++k) {
		getline(file_in,kline); // read k-th line
		vector<string> x_str=split(kline,' '); // split line
		vector<int> neighbs_k;
		int cnt=0;
		for (size_t kk=0;kk<x_str.size();++kk) {
			if (not x_str[kk].empty()) {
				cnt+=1;
				if (cnt>3) {
					neighbs_k.push_back(atoi(x_str[kk].c_str())-1);
				}
			}
		}
		neighbs.push_back(neighbs_k);
	}
	//
	file_in.close();
	//
	return neighbs; 
}


//////////////////////////////// create cliques with metadata for bnd_resolve 
// data structure for a marked point
struct marked_point {
	vector<double> x;
	double a,b,th;
	vector<double> na,nb,z,zsqrt;
};
// create a marked point
marked_point create_mp(double th, double a, double b, vector<double> x) {
	vector<double> na(2),nb(2),z(3),zsqrt(3);
	na[0]=cos(th);
	na[1]=sin(th);	
	nb[0]=sin(th);
	nb[1]=-cos(th);
	z[0]=pow((na[0]/a),2)+pow((nb[0]/b),2);
	z[1]=pow((na[1]/a),2)+pow((nb[1]/b),2);
	z[2]=na[0]*na[1]/pow(a,2)+nb[0]*nb[1]/pow(b,2);
	zsqrt[0]=na[0]*na[0]/a+nb[0]*nb[0]/b;
	zsqrt[1]=na[1]*na[1]/a+nb[1]*nb[1]/b;
	zsqrt[2]=na[0]*na[1]/a+nb[0]*nb[1]/b;
	// change nucleation pt to handle periodicity (if necessary)
	return {x,a,b,th,na,nb,z,zsqrt};
}
// data structure for a clique
struct clique {
	vector<marked_point> mp;
	vector<int> neighbs;
	//
	int igrain_id;
	vector<int> lorder;
	vector<double> c_ij;
	vector<vector<vector<vector<double>>>> thsr_ij;
	vector<double> cmn_thsr;
	vector<vector<int>> cmn_types;
	vector<vector<double>> uij1,uij2,zij,xij;	
};
// create list of all cliques solved by bnd_solve 
vector<clique> create_cliques(bnd_solve_data data,marked_point_pattern mpp,vector<vector<int>> neighbs) {
	vector<clique> cliques;
	for (int i=0;i<data.ngrains;++i) {
		int k=data.igrain_ids[i];
		vector<marked_point> mp_vec={create_mp(mpp.th[k],mpp.a[k],mpp.b[k],mpp.x[k])};
		vector<double> x0=mpp.x[k];
		//
		// handle periodicity by changing nucleation points of neighbors in clique 
		for (size_t j=0;j<neighbs[k].size();++j) {
			vector<double> xj=mpp.x[neighbs[k][j]];
			if ((xj[0]-x0[0])>mpp.lx/2.) {
				xj[0]-=mpp.lx;
			}
			else if ((xj[0]-x0[0])<-mpp.lx/2.) {
				xj[0]+=mpp.lx;
			}
			if ((xj[1]-x0[1])>mpp.ly/2.) {
				xj[1]-=mpp.ly;
			}
			else if ((xj[1]-x0[1])<-mpp.ly/2.) {
				xj[1]+=mpp.ly;
			}
			mp_vec.push_back(create_mp(mpp.th[neighbs[k][j]],mpp.a[neighbs[k][j]],mpp.b[neighbs[k][j]],xj));
		}
		clique clique_i;
		clique_i.mp=mp_vec;
		clique_i.neighbs=neighbs[k];
		clique_i.igrain_id=data.igrain_ids[i];
		clique_i.lorder=data.lorder[i];
		clique_i.c_ij=data.c_ij[i];
		clique_i.thsr_ij=data.thsr_ij[i];	
		clique_i.cmn_thsr=data.cmn_thsr[i];
		clique_i.cmn_types=data.cmn_types[i];
		clique_i.uij1=data.uij1[i];
		clique_i.uij2=data.uij2[i];
		clique_i.zij=data.zij[i];
		clique_i.xij=data.xij[i];	
		//
		cliques.push_back(clique_i);
	}
	return cliques;
}

////////////////////////////////  core of bnd_resolve 
//
// first parameterization of S^1
vector<double> xth1(double th,vector<double> uij1,vector<double> uij2) {
	return {uij1[0]*cos(th)+uij2[0]*sin(th),uij1[1]*cos(th)+uij2[1]*sin(th)};
}
// Corresponding derivative
vector<double> dxth1(double th,vector<double> uij1,vector<double> uij2) {
	return {-uij1[0]*sin(th)+uij2[0]*cos(th),-uij1[1]*sin(th)+uij2[1]*cos(th)};
}
// second parameterization of S^1
vector<double> xth2(double th,vector<double> xij) {
	vector<vector<double>> rth;
	rth.push_back({cos(th),-sin(th)});
	rth.push_back({sin(th),cos(th)});	
	vector<double> x(2);
	x[0]=(rth[0][0]*xij[0]+rth[0][1]*xij[1])/sqrt(xij[0]*xij[0]+xij[1]*xij[1]);
	x[1]=(rth[1][0]*xij[0]+rth[1][1]*xij[1])/sqrt(xij[0]*xij[0]+xij[1]*xij[1]);
	return x;
}
// Corresponding derivative
vector<double> dxth2(double th,vector<double> xij) {
	vector<vector<double>> drth;
	drth.push_back({-sin(th),-cos(th)});
	drth.push_back({cos(th),-sin(th)});	
	vector<double> dx(2);
	dx[0]=(drth[0][0]*xij[0]+drth[0][1]*xij[1])/sqrt(xij[0]*xij[0]+xij[1]*xij[1]);
	dx[1]=(drth[1][0]*xij[0]+drth[1][1]*xij[1])/sqrt(xij[0]*xij[0]+xij[1]*xij[1]);
	return dx;	
}
// data structure of parameterization for a cmn curve
struct cmn_curve_param {
	clique clq;
	int t;
	vector<double> na,nb,x,zsqrt,zsqrt_inv;
	double om_t;
	vector<double> zij,xij;
	//	
	void set_clq(clique _clq) {
		clq=_clq;
		na=clq.mp[0].na;
		nb=clq.mp[0].nb;
		x=clq.mp[0].x;
		zsqrt=clq.mp[0].zsqrt;
		double det=zsqrt[0]*zsqrt[1]-zsqrt[2]*zsqrt[2];
		zsqrt_inv={zsqrt[1]/det,zsqrt[0]/det,-zsqrt[2]/det};
	}
	//
	void set_t(int _t) {
		t=_t;
		if (clq.c_ij[t]>=1) {
			om_t=sign(clq.uij1[t][0]*na[0]+clq.uij1[t][1]*na[1])*acos(clq.uij1[t][0]*nb[0]+clq.uij1[t][1]*nb[1]);
			}
		else {
			om_t=sign(na[0]*clq.xij[t][0]+na[1]*clq.xij[t][1])*acos((nb[0]*clq.xij[t][0]+nb[1]*clq.xij[t][1])/sqrt(clq.xij[t][0]*clq.xij[t][0]+clq.xij[t][1]*clq.xij[t][1]));
		}
		zij=clq.zij[t];
		xij=clq.xij[t];
	}
	// global parameterization of S^1
	vector<double> xth(double th) {
		if (clq.c_ij[t]>=1) {
			return xth1(th-om_t,clq.uij1[t],clq.uij2[t]);
		}
		else {
			return xth2(th-om_t,clq.xij[t]);
		}	
	}
	// first derivative of global parameterization of S^1
	vector<double> dxth(double th) {
		if (clq.c_ij[t]>=1) {
			return dxth1(th-om_t,clq.uij1[t],clq.uij2[t]);
		}
		else {
			return dxth2(th-om_t,clq.xij[t]);
		}
	}
	// contact function (1st type)
	double xi_t0(double th) {
		vector<double> x=xth(th);
		double a=xij[0]*xij[0]*zij[0]+xij[1]*xij[1]*zij[1]+2.*xij[0]*xij[1]*zij[2];
		double b=x[0]*xij[0]*zij[0]+x[1]*xij[1]*zij[1]+(x[0]*xij[1]+x[1]*xij[0])*zij[2];
		double c=x[0]*x[0]*zij[0]+x[1]*x[1]*zij[1]+2.*x[0]*x[1]*zij[2];
		//
		return a/(b+sqrt(b*b-a*(c-1.)));
	}
	// derivative of contact function (1st type)
	double dxi_t0(double th) {
		vector<double> x=xth(th);
		vector<double> dx=dxth(th);
		double a=xij[0]*xij[0]*zij[0]+xij[1]*xij[1]*zij[1]+2.*xij[0]*xij[1]*zij[2];
		double b=x[0]*xij[0]*zij[0]+x[1]*xij[1]*zij[1]+(x[0]*xij[1]+x[1]*xij[0])*zij[2];
		double c=x[0]*x[0]*zij[0]+x[1]*x[1]*zij[1]+2.*x[0]*x[1]*zij[2];
		double d=dx[0]*xij[0]*zij[0]+dx[1]*xij[1]*zij[1]+(dx[0]*xij[1]+dx[1]*xij[0])*zij[2];
		double e=dx[0]*x[0]*zij[0]+dx[1]*x[1]*zij[1]+(dx[0]*x[1]+dx[1]*x[0])*zij[2];
		//
		return -pow(xi_t0(th),2)/a*(d+(b*d-a*e)/sqrt(b*b-a*(c-1.)));
	}
	// contact function (2nd type)
	double xi_t1(double th) {
		vector<double> x=xth(th);
		double a=xij[0]*xij[0]*zij[0]+xij[1]*xij[1]*zij[1]+2.*xij[0]*xij[1]*zij[2];
		double b=x[0]*xij[0]*zij[0]+x[1]*xij[1]*zij[1]+(x[0]*xij[1]+x[1]*xij[0])*zij[2];
		double c=x[0]*x[0]*zij[0]+x[1]*x[1]*zij[1]+2.*x[0]*x[1]*zij[2];
		//
		return a/(b-sqrt(b*b-a*(c-1.)));
	}
	// derivative of contact function (2nd type)
	double dxi_t1(double th) {
		vector<double> x=xth(th);
		vector<double> dx=dxth(th);
		double a=xij[0]*xij[0]*zij[0]+xij[1]*xij[1]*zij[1]+2.*xij[0]*xij[1]*zij[2];
		double b=x[0]*xij[0]*zij[0]+x[1]*xij[1]*zij[1]+(x[0]*xij[1]+x[1]*xij[0])*zij[2];
		double c=x[0]*x[0]*zij[0]+x[1]*x[1]*zij[1]+2.*x[0]*x[1]*zij[2];
		double d=dx[0]*xij[0]*zij[0]+dx[1]*xij[1]*zij[1]+(dx[0]*xij[1]+dx[1]*xij[0])*zij[2];
		double e=dx[0]*x[0]*zij[0]+dx[1]*x[1]*zij[1]+(dx[0]*x[1]+dx[1]*x[0])*zij[2];
		//
		return -pow(xi_t0(th),2)/a*(d-(b*d-a*e)/sqrt(b*b-a*(c-1.)));
	}
	// 
	double x_t0(double th) {
		vector<vector<double>> rth;
		rth.push_back({cos(th),-sin(th)});
		rth.push_back({sin(th),cos(th)});
		return x[0]+xi_t0(th)*(zsqrt_inv[0]*(rth[0][0]*nb[0]+rth[0][1]*nb[1])+zsqrt_inv[2]*(rth[1][0]*nb[0]+rth[1][1]*nb[1]));
	}
	// 
	double y_t0(double th) {
		vector<vector<double>> rth;
		rth.push_back({cos(th),-sin(th)});
		rth.push_back({sin(th),cos(th)});	
		return x[1]+xi_t0(th)*(zsqrt_inv[2]*(rth[0][0]*nb[0]+rth[0][1]*nb[1])+zsqrt_inv[1]*(rth[1][0]*nb[0]+rth[1][1]*nb[1]));
	}
	// 
	double x_t1(double th) {
		vector<vector<double>> rth;
		rth.push_back({cos(th),-sin(th)});
		rth.push_back({sin(th),cos(th)});	
		return x[0]+xi_t1(th)*(zsqrt_inv[0]*(rth[0][0]*nb[0]+rth[0][1]*nb[1])+zsqrt_inv[2]*(rth[1][0]*nb[0]+rth[1][1]*nb[1]));
	}
	// 
	double y_t1(double th) {
		vector<vector<double>> rth;
		rth.push_back({cos(th),-sin(th)});
		rth.push_back({sin(th),cos(th)});	
		return x[1]+xi_t1(th)*(zsqrt_inv[2]*(rth[0][0]*nb[0]+rth[0][1]*nb[1])+zsqrt_inv[1]*(rth[1][0]*nb[0]+rth[1][1]*nb[1]));
	}	
};


// data structure for solved bnd
struct bnd {
	vector<vector<double>> x,y,th;
	vector<vector<double>> x_w1,y_w1,x_w0;
	vector<vector<double>> x_w2,y_w2,y_w0;
	vector<vector<double>> xi,dxi;
	vector<vector<double>> L,nx,ny,dt;
	double W0;
	vector<double> W0_10,W0_20,W0_30,W0_40,W0_50,W0_60;
	double W1;
	vector<double> W1_10,W1_20,W1_30,W1_40,W1_50,W1_60;
	vector<double> W1_01,W1_02,W1_03,W1_04,W1_05,W1_06;
	vector<double> W1_11;
	vector<double> W1_21,W1_12;
	vector<double> W1_31,W1_22,W1_13;
	vector<double> W1_41,W1_32,W1_23,W1_14;
	vector<double> W1_51,W1_42,W1_33,W1_24,W1_15;
	double W2;
	vector<double> W2_10,W2_20,W2_30,W2_40,W2_50,W2_60;
	vector<double> W2_01,W2_02,W2_03,W2_04,W2_05,W2_06;
	vector<double> W2_11;
	vector<double> W2_21,W2_12;
	vector<double> W2_31,W2_22,W2_13;
	vector<double> W2_41,W2_32,W2_23,W2_14;
	vector<double> W2_51,W2_42,W2_33,W2_24,W2_15;

};
vector<bnd> resolve_bnds(vector<clique> cliques, int npts) {
	vector<bnd> bnds;
	//
	for (size_t i=0;i<cliques.size();++i) {
		cout << "working on clique # " << i << endl;
		clique clq=cliques[i];
		bnd bnd_i;
		bnd_i.W0=0.;
		bnd_i.W0_10={0.,0.};
		bnd_i.W0_20={0.,0.,0.};
		bnd_i.W0_30={0.,0.,0.,0.};
		bnd_i.W0_40={0.,0.,0.,0.,0.};
		bnd_i.W0_50={0.,0.,0.,0.,0.,0.};
		bnd_i.W0_60={0.,0.,0.,0.,0.,0.,0.};
		bnd_i.W1=0.;
		bnd_i.W1_10={0.,0.};bnd_i.W1_01={0.,0.};
		bnd_i.W1_20={0.,0.,0.};bnd_i.W1_11={0.,0.,0.};bnd_i.W1_02={0.,0.,0.};
		bnd_i.W1_30={0.,0.,0.,0.};bnd_i.W1_21={0.,0.,0.,0.};bnd_i.W1_12={0.,0.,0.,0.};bnd_i.W1_03={0.,0.,0.,0.};
		bnd_i.W1_40={0.,0.,0.,0.,0.};bnd_i.W1_31={0.,0.,0.,0.,0.};bnd_i.W1_22={0.,0.,0.,0.,0.};bnd_i.W1_13={0.,0.,0.,0.,0.};bnd_i.W1_04={0.,0.,0.,0.,0.};
		bnd_i.W1_50={0.,0.,0.,0.,0.,0.};bnd_i.W1_41={0.,0.,0.,0.,0.,0.};bnd_i.W1_32={0.,0.,0.,0.,0.,0.};bnd_i.W1_23={0.,0.,0.,0.,0.,0.};bnd_i.W1_14={0.,0.,0.,0.,0.,0.};bnd_i.W1_05={0.,0.,0.,0.,0.,0.};
		bnd_i.W1_60={0.,0.,0.,0.,0.,0.,0.};bnd_i.W1_51={0.,0.,0.,0.,0.,0.,0.};bnd_i.W1_42={0.,0.,0.,0.,0.,0.,0.};bnd_i.W1_33={0.,0.,0.,0.,0.,0.,0.};bnd_i.W1_24={0.,0.,0.,0.,0.,0.,0.};bnd_i.W1_15={0.,0.,0.,0.,0.,0.,0.};bnd_i.W1_06={0.,0.,0.,0.,0.,0.,0.};
		bnd_i.W2=0.;
		bnd_i.W2_10={0.,0.};bnd_i.W2_01={0.,0.};
		bnd_i.W2_20={0.,0.,0.};bnd_i.W2_11={0.,0.,0.};bnd_i.W2_02={0.,0.,0.};
		bnd_i.W2_30={0.,0.,0.,0.};bnd_i.W2_21={0.,0.,0.,0.};bnd_i.W2_12={0.,0.,0.,0.};bnd_i.W2_03={0.,0.,0.,0.};
		bnd_i.W2_40={0.,0.,0.,0.,0.};bnd_i.W2_31={0.,0.,0.,0.,0.};bnd_i.W2_22={0.,0.,0.,0.,0.};bnd_i.W2_13={0.,0.,0.,0.,0.};bnd_i.W2_04={0.,0.,0.,0.,0.};
		bnd_i.W2_50={0.,0.,0.,0.,0.,0.};bnd_i.W2_41={0.,0.,0.,0.,0.,0.};bnd_i.W2_32={0.,0.,0.,0.,0.,0.};bnd_i.W2_23={0.,0.,0.,0.,0.,0.};bnd_i.W2_14={0.,0.,0.,0.,0.,0.};bnd_i.W2_05={0.,0.,0.,0.,0.,0.};
		bnd_i.W2_60={0.,0.,0.,0.,0.,0.,0.};bnd_i.W2_51={0.,0.,0.,0.,0.,0.,0.};bnd_i.W2_42={0.,0.,0.,0.,0.,0.,0.};bnd_i.W2_33={0.,0.,0.,0.,0.,0.,0.};bnd_i.W2_24={0.,0.,0.,0.,0.,0.,0.};bnd_i.W2_15={0.,0.,0.,0.,0.,0.,0.};bnd_i.W2_06={0.,0.,0.,0.,0.,0.,0.};
		//
		vector<double> na=clq.mp[0].na;
		vector<double> nb=clq.mp[0].nb;
		vector<double> x=clq.mp[0].x;
		vector<double> zsqrt=clq.mp[0].zsqrt;
		//
		//int igrain_id=clq.igrain_id;
		vector<int> lorder=clq.lorder;
		vector<double> c_ij=clq.c_ij;
		vector<vector<vector<vector<double>>>> thsr_ij=clq.thsr_ij;
		vector<double> cmn_thsr=clq.cmn_thsr;
		vector<vector<int>> cmn_types=clq.cmn_types;
		vector<vector<double>> uij1=clq.uij1;
		vector<vector<double>> uij2=clq.uij2;
		vector<vector<double>> zij=clq.zij;
		vector<vector<double>> xij=clq.xij;
		
		
		
		vector<double> bnd_i_W0_10(2);
		
		
		//
		for (size_t j=0;j<lorder.size()-1;++j) {
			int t=lorder[j];
			int s=cmn_types.size();
			cmn_curve_param prm;
			prm.set_clq(clq);
			prm.set_t(t);
			// had to create re_j(j,s) because cmn_types[-1] (and cmn_thsr[-1]) is not a valid c++ 
			// indexation for the last element of cmn_types (and cmn_thsr).
			//
			// first case: boundary completely radially convex
			

			if ((cmn_types[re_j(j,s)-1][1]==0)&(cmn_types[j][0]==0)) {
//				cout << "test 1a" << endl;
				int loc_cnt=1;
				if (thsr_ij[lorder[j]][0].size()==2) {
					if ((cmn_thsr[re_j(j,s)-1]<thsr_ij[lorder[j]][0][0][1])&(cmn_thsr[re_j(j,s)-1]>thsr_ij[lorder[j]][0][0][0])&(cmn_thsr[j]<thsr_ij[lorder[j]][0][0][1])&(cmn_thsr[j]>thsr_ij[lorder[j]][0][0][0])) {
						bnd_i.th.push_back(linspace(cmn_thsr[re_j(j,s)-1],cmn_thsr[j],npts));
					}
					else if ((cmn_thsr[re_j(j,s)-1]<thsr_ij[lorder[j]][0][1][1])&(cmn_thsr[re_j(j,s)-1]>thsr_ij[lorder[j]][0][1][0])&(cmn_thsr[j]<thsr_ij[lorder[j]][0][1][1])&(cmn_thsr[j]>thsr_ij[lorder[j]][0][1][0])) {
						bnd_i.th.push_back(linspace(cmn_thsr[re_j(j,s)-1],cmn_thsr[j],npts));
					}
					else {
						loc_cnt=2;
						bnd_i.th.push_back(linspace(cmn_thsr[re_j(j,s)-1],thsr_ij[lorder[j]][0][1][1],npts));
						bnd_i.th.push_back(linspace(thsr_ij[lorder[j]][0][0][0],cmn_thsr[j],npts));
					}
				}
				else {
					bnd_i.th.push_back(linspace(cmn_thsr[re_j(j,s)-1],cmn_thsr[j],npts));
				}
//				cout << "test 1b" << endl;
				for (int i_cnt=loc_cnt;i_cnt>0;--i_cnt) {
					int ss=bnd_i.th.size();
					int n=bnd_i.th[ss-i_cnt].size();
					vector<double> vec_xi(n),vec_dxi(n),vec_x(n),vec_y(n);
					vector<double> vec_nx(n-1),vec_ny(n-1),vec_L(n-1),vec_dt(n-1);
					for (int k=0;k<n;++k) {
						vec_xi[k]=prm.xi_t0(bnd_i.th[ss-i_cnt][k]);
						vec_dxi[k]=prm.dxi_t0(bnd_i.th[ss-i_cnt][k]);
						vec_x[k]=prm.x_t0(bnd_i.th[ss-i_cnt][k]);
						vec_y[k]=prm.y_t0(bnd_i.th[ss-i_cnt][k]);
//						cout << "test 1d" << endl;
						if (k>0) {
							vec_L[k-1]=pow(pow(vec_x[k]-vec_x[k-1],2.)+pow(vec_y[k]-vec_y[k-1],2.),.5);
							bnd_i.W1+=vec_L[k-1]/2.;
							vec_nx[k-1]=-(vec_y[k-1]-vec_y[k])/vec_L[k-1];
							vec_ny[k-1]=(vec_x[k-1]-vec_x[k])/vec_L[k-1];
							W_rs_nu dW1_10=get_W1_rs_type_II_k(1,0,{vec_x[k-1],vec_y[k-1]},{vec_x[k],vec_y[k]},{vec_nx[k-1],vec_ny[k-1]},vec_L[k-1]);
							bnd_i.W1_10[0]+=dW1_10.comps[0];
							bnd_i.W1_10[1]+=dW1_10.comps[1];
							W_rs_nu dW1_01=get_W1_rs_type_II_k(0,1,{vec_x[k-1],vec_y[k-1]},{vec_x[k],vec_y[k]},{vec_nx[k-1],vec_ny[k-1]},vec_L[k-1]);
							bnd_i.W1_01[0]+=dW1_01.comps[0];
							bnd_i.W1_01[1]+=dW1_01.comps[1];
							bnd_i.W0+=1./4.*vec_L[k-1]*((vec_x[k-1]+vec_x[k])*vec_nx[k-1]+(vec_y[k-1]+vec_y[k])*vec_ny[k-1]);	
							W_rs_nu dW0_10=get_W0_rs_type_II_k(1,0,{vec_x[k-1],vec_y[k-1]},{vec_x[k],vec_y[k]},{vec_nx[k-1],vec_ny[k-1]},vec_L[k-1]);
							bnd_i.W0_10[0]+=dW0_10.comps[0];
							bnd_i.W0_10[1]+=dW0_10.comps[1];
										
						}
						if ((k==1)&(bnd_i.L.size()>0)) {
							vec_dt[0]=sign(bnd_i.nx.back().back()*vec_ny[0]-bnd_i.ny.back().back()*vec_nx[0])*acos(bnd_i.nx.back().back()*vec_nx[0]+bnd_i.ny.back().back()*vec_ny[0]);
							bnd_i.W2+=.5*int_k(0,0,vec_dt[0]);
							W_rs_nu dW2_10=get_W2_rs_type_II_k(1,0,{bnd_i.x.back().end()[-2],bnd_i.y.back().end()[-2]},{bnd_i.x.back().back(),bnd_i.y.back().back()},{bnd_i.nx.back().back(),bnd_i.ny.back().back()},bnd_i.L.back().back(),vec_dt[k-1]);
							bnd_i.W2_10[0]+=dW2_10.comps[0];
							bnd_i.W2_10[1]+=dW2_10.comps[1];	
							W_rs_nu dW2_01=get_W2_rs_type_II_k(0,1,{bnd_i.x.back().end()[-2],bnd_i.y.back().end()[-2]},{bnd_i.x.back().back(),bnd_i.y.back().back()},{bnd_i.nx.back().back(),bnd_i.ny.back().back()},bnd_i.L.back().back(),vec_dt[k-1]);
							bnd_i.W2_01[0]+=dW2_01.comps[0];
							bnd_i.W2_01[1]+=dW2_01.comps[1];		
						} 
						else if (k>1) {
							vec_dt[k-1]=sign(vec_nx[k-2]*vec_ny[k-1]-vec_ny[k-2]*vec_nx[k-1])*acos(vec_nx[k-2]*vec_nx[k-1]+vec_ny[k-2]*vec_ny[k-1]);
							bnd_i.W2+=.5*int_k(0,0,vec_dt[k-1]);
							W_rs_nu dW2_10=get_W2_rs_type_II_k(1,0,{vec_x[k-2],vec_y[k-2]},{vec_x[k-1],vec_y[k-1]},{vec_nx[k-2],vec_ny[k-2]},vec_L[k-2],vec_dt[k-1]);
							bnd_i.W2_10[0]+=dW2_10.comps[0];
							bnd_i.W2_10[1]+=dW2_10.comps[1];
							W_rs_nu dW2_01=get_W2_rs_type_II_k(0,1,{vec_x[k-2],vec_y[k-2]},{vec_x[k-1],vec_y[k-1]},{vec_nx[k-2],vec_ny[k-2]},vec_L[k-2],vec_dt[k-1]);
							bnd_i.W2_01[0]+=dW2_01.comps[0];
							bnd_i.W2_01[1]+=dW2_01.comps[1];
						}
					}
					bnd_i.xi.push_back(vec_xi);
					bnd_i.dxi.push_back(vec_dxi);
					bnd_i.x.push_back(vec_x);
					bnd_i.y.push_back(vec_y);
//					cout << "test 1e" << endl;
					bnd_i.L.push_back(vec_L);
					bnd_i.nx.push_back(vec_nx);
					bnd_i.ny.push_back(vec_ny);
					bnd_i.dt.push_back(vec_dt);
				}
//				cout << "test 1f" << endl;
			}
			//
			// second case
			else if ((cmn_types[re_j(j,s)-1][1]==1)&(cmn_types[j][0]==0)) {
//				cout << "test 2a" << endl;
				int loc_cnt=1;
				if (thsr_ij[lorder[j]][1].size()==2) {
					if ((cmn_thsr[re_j(j,s)-1]<thsr_ij[lorder[j]][1][0][1])&(cmn_thsr[re_j(j,s)-1]>thsr_ij[lorder[j]][1][0][0])) {
						loc_cnt=2;
						bnd_i.th.push_back(linspace(cmn_thsr[re_j(j,s)-1],thsr_ij[lorder[j]][1][0][0],npts));
						bnd_i.th.push_back(linspace(thsr_ij[lorder[j]][1][1][1],thsr_ij[lorder[j]][1][1][0],npts));
					}
					else if ((cmn_thsr[re_j(j,s)-1]<thsr_ij[lorder[j]][1][1][1])&(cmn_thsr[re_j(j,s)-1]>thsr_ij[lorder[j]][1][1][0])) {
						bnd_i.th.push_back(linspace(cmn_thsr[re_j(j,s)-1],thsr_ij[lorder[j]][1][1][0],npts));
					}
				}
				else {
					if ((cmn_thsr[re_j(j,s)-1]<thsr_ij[lorder[j]][1][0][1])&(cmn_thsr[re_j(j,s)-1]>thsr_ij[lorder[j]][1][0][0])) {
						bnd_i.th.push_back(linspace(cmn_thsr[re_j(j,s)-1],thsr_ij[lorder[j]][1][0][0],npts));
					}
				}
				for (int i_cnt=loc_cnt;i_cnt>0;--i_cnt) {
					int ss=bnd_i.th.size();
					int n=bnd_i.th[ss-i_cnt].size();
					vector<double> vec_xi(n),vec_dxi(n),vec_x(n),vec_y(n);
					vector<double> vec_nx(n-1),vec_ny(n-1),vec_L(n-1),vec_dt(n-1);
					for (int k=0;k<n;++k) {
						vec_xi[k]=prm.xi_t1(bnd_i.th[ss-i_cnt][k]);
						vec_dxi[k]=prm.dxi_t1(bnd_i.th[ss-i_cnt][k]);
						vec_x[k]=prm.x_t1(bnd_i.th[ss-i_cnt][k]);
						vec_y[k]=prm.y_t1(bnd_i.th[ss-i_cnt][k]);
						//
						if (k>0) {
							vec_L[k-1]=pow(pow(vec_x[k]-vec_x[k-1],2.)+pow(vec_y[k]-vec_y[k-1],2.),.5);
							bnd_i.W1+=vec_L[k-1]/2.;
							vec_nx[k-1]=-(vec_y[k-1]-vec_y[k])/vec_L[k-1];
							vec_ny[k-1]=(vec_x[k-1]-vec_x[k])/vec_L[k-1];
							W_rs_nu dW1_10=get_W1_rs_type_II_k(1,0,{vec_x[k-1],vec_y[k-1]},{vec_x[k],vec_y[k]},{vec_nx[k-1],vec_ny[k-1]},vec_L[k-1]);
							bnd_i.W1_10[0]+=dW1_10.comps[0];
							bnd_i.W1_10[1]+=dW1_10.comps[1];
							W_rs_nu dW1_01=get_W1_rs_type_II_k(0,1,{vec_x[k-1],vec_y[k-1]},{vec_x[k],vec_y[k]},{vec_nx[k-1],vec_ny[k-1]},vec_L[k-1]);
							bnd_i.W1_01[0]+=dW1_01.comps[0];
							bnd_i.W1_01[1]+=dW1_01.comps[1];
							bnd_i.W0+=1./4.*vec_L[k-1]*((vec_x[k-1]+vec_x[k])*vec_nx[k-1]+(vec_y[k-1]+vec_y[k])*vec_ny[k-1]);	
							W_rs_nu dW0_10=get_W0_rs_type_II_k(1,0,{vec_x[k-1],vec_y[k-1]},{vec_x[k],vec_y[k]},{vec_nx[k-1],vec_ny[k-1]},vec_L[k-1]);
							bnd_i.W0_10[0]+=dW0_10.comps[0];
							bnd_i.W0_10[1]+=dW0_10.comps[1];
						}	
						if ((k==1)&(bnd_i.L.size()>0)) {
							vec_dt[0]=sign(bnd_i.nx.back().back()*vec_ny[0]-bnd_i.ny.back().back()*vec_nx[0])*acos(bnd_i.nx.back().back()*vec_nx[0]+bnd_i.ny.back().back()*vec_ny[0]);
							bnd_i.W2+=.5*int_k(0,0,vec_dt[0]);
							W_rs_nu dW2_10=get_W2_rs_type_II_k(1,0,{bnd_i.x.back().end()[-2],bnd_i.y.back().end()[-2]},{bnd_i.x.back().back(),bnd_i.y.back().back()},{bnd_i.nx.back().back(),bnd_i.ny.back().back()},bnd_i.L.back().back(),vec_dt[k-1]);
							bnd_i.W2_10[0]+=dW2_10.comps[0];
							bnd_i.W2_10[1]+=dW2_10.comps[1];	
							W_rs_nu dW2_01=get_W2_rs_type_II_k(0,1,{bnd_i.x.back().end()[-2],bnd_i.y.back().end()[-2]},{bnd_i.x.back().back(),bnd_i.y.back().back()},{bnd_i.nx.back().back(),bnd_i.ny.back().back()},bnd_i.L.back().back(),vec_dt[k-1]);
							bnd_i.W2_01[0]+=dW2_01.comps[0];
							bnd_i.W2_01[1]+=dW2_01.comps[1];		
						} 
						else if (k>1) {
							vec_dt[k-1]=sign(vec_nx[k-2]*vec_ny[k-1]-vec_ny[k-2]*vec_nx[k-1])*acos(vec_nx[k-2]*vec_nx[k-1]+vec_ny[k-2]*vec_ny[k-1]);
							bnd_i.W2+=.5*int_k(0,0,vec_dt[k-1]);
							W_rs_nu dW2_10=get_W2_rs_type_II_k(1,0,{vec_x[k-2],vec_y[k-2]},{vec_x[k-1],vec_y[k-1]},{vec_nx[k-2],vec_ny[k-2]},vec_L[k-2],vec_dt[k-1]);
							bnd_i.W2_10[0]+=dW2_10.comps[0];
							bnd_i.W2_10[1]+=dW2_10.comps[1];
							W_rs_nu dW2_01=get_W2_rs_type_II_k(0,1,{vec_x[k-2],vec_y[k-2]},{vec_x[k-1],vec_y[k-1]},{vec_nx[k-2],vec_ny[k-2]},vec_L[k-2],vec_dt[k-1]);
							bnd_i.W2_01[0]+=dW2_01.comps[0];
							bnd_i.W2_01[1]+=dW2_01.comps[1];
						}
					}
					bnd_i.xi.push_back(vec_xi);
					bnd_i.dxi.push_back(vec_dxi);
					bnd_i.x.push_back(vec_x);
					bnd_i.y.push_back(vec_y);
					//
					bnd_i.L.push_back(vec_L);
					bnd_i.nx.push_back(vec_nx);
					bnd_i.ny.push_back(vec_ny);
					bnd_i.dt.push_back(vec_dt);
				}				
				//
				loc_cnt=1;
				if (thsr_ij[lorder[j]][0].size()==2) {
					if ((cmn_thsr[j]<thsr_ij[lorder[j]][0][0][1])&(cmn_thsr[j]>thsr_ij[lorder[j]][0][0][0])) {
						loc_cnt=2;
						bnd_i.th.push_back(linspace(thsr_ij[lorder[j]][0][1][0],thsr_ij[lorder[j]][0][1][1],npts));
						bnd_i.th.push_back(linspace(thsr_ij[lorder[j]][0][0][0],cmn_thsr[j],npts));
					}
					else if ((cmn_thsr[j]<thsr_ij[lorder[j]][0][1][1])&(cmn_thsr[j]>thsr_ij[lorder[j]][0][1][0])) {
						bnd_i.th.push_back(linspace(thsr_ij[lorder[j]][0][1][0],cmn_thsr[j],npts));
					}
				}
				else {
					if ((cmn_thsr[j]<thsr_ij[lorder[j]][0][0][1])&(cmn_thsr[j]>thsr_ij[lorder[j]][0][0][0])) {
						bnd_i.th.push_back(linspace(thsr_ij[lorder[j]][0][0][0],cmn_thsr[j],npts));
					}
				}
				for (int i_cnt=loc_cnt;i_cnt>0;--i_cnt) {
					int ss=bnd_i.th.size();
					int n=bnd_i.th[ss-i_cnt].size();
					vector<double> vec_xi(n),vec_dxi(n),vec_x(n),vec_y(n);
					vector<double> vec_nx(n-1),vec_ny(n-1),vec_L(n-1),vec_dt(n-1);
					for (int k=0;k<n;++k) {
						vec_xi[k]=prm.xi_t0(bnd_i.th[ss-i_cnt][k]);
						vec_dxi[k]=prm.dxi_t0(bnd_i.th[ss-i_cnt][k]);
						vec_x[k]=prm.x_t0(bnd_i.th[ss-i_cnt][k]);
						vec_y[k]=prm.y_t0(bnd_i.th[ss-i_cnt][k]);
						//
						if (k>0) {
							vec_L[k-1]=pow(pow(vec_x[k]-vec_x[k-1],2.)+pow(vec_y[k]-vec_y[k-1],2.),.5);
							bnd_i.W1+=vec_L[k-1]/2.;
							vec_nx[k-1]=-(vec_y[k-1]-vec_y[k])/vec_L[k-1];
							vec_ny[k-1]=(vec_x[k-1]-vec_x[k])/vec_L[k-1];
							W_rs_nu dW1_10=get_W1_rs_type_II_k(1,0,{vec_x[k-1],vec_y[k-1]},{vec_x[k],vec_y[k]},{vec_nx[k-1],vec_ny[k-1]},vec_L[k-1]);
							bnd_i.W1_10[0]+=dW1_10.comps[0];
							bnd_i.W1_10[1]+=dW1_10.comps[1];
							W_rs_nu dW1_01=get_W1_rs_type_II_k(0,1,{vec_x[k-1],vec_y[k-1]},{vec_x[k],vec_y[k]},{vec_nx[k-1],vec_ny[k-1]},vec_L[k-1]);
							bnd_i.W1_01[0]+=dW1_01.comps[0];
							bnd_i.W1_01[1]+=dW1_01.comps[1];
							bnd_i.W0+=1./4.*vec_L[k-1]*((vec_x[k-1]+vec_x[k])*vec_nx[k-1]+(vec_y[k-1]+vec_y[k])*vec_ny[k-1]);	
							W_rs_nu dW0_10=get_W0_rs_type_II_k(1,0,{vec_x[k-1],vec_y[k-1]},{vec_x[k],vec_y[k]},{vec_nx[k-1],vec_ny[k-1]},vec_L[k-1]);
							bnd_i.W0_10[0]+=dW0_10.comps[0];
							bnd_i.W0_10[1]+=dW0_10.comps[1];
						}	
						if ((k==1)&(bnd_i.L.size()>0)) {
							vec_dt[0]=sign(bnd_i.nx.back().back()*vec_ny[0]-bnd_i.ny.back().back()*vec_nx[0])*acos(bnd_i.nx.back().back()*vec_nx[0]+bnd_i.ny.back().back()*vec_ny[0]);
							bnd_i.W2+=.5*int_k(0,0,vec_dt[0]);
							W_rs_nu dW2_10=get_W2_rs_type_II_k(1,0,{bnd_i.x.back().end()[-2],bnd_i.y.back().end()[-2]},{bnd_i.x.back().back(),bnd_i.y.back().back()},{bnd_i.nx.back().back(),bnd_i.ny.back().back()},bnd_i.L.back().back(),vec_dt[k-1]);
							bnd_i.W2_10[0]+=dW2_10.comps[0];
							bnd_i.W2_10[1]+=dW2_10.comps[1];	
							W_rs_nu dW2_01=get_W2_rs_type_II_k(0,1,{bnd_i.x.back().end()[-2],bnd_i.y.back().end()[-2]},{bnd_i.x.back().back(),bnd_i.y.back().back()},{bnd_i.nx.back().back(),bnd_i.ny.back().back()},bnd_i.L.back().back(),vec_dt[k-1]);
							bnd_i.W2_01[0]+=dW2_01.comps[0];
							bnd_i.W2_01[1]+=dW2_01.comps[1];		
						} 
						else if (k>1) {
							vec_dt[k-1]=sign(vec_nx[k-2]*vec_ny[k-1]-vec_ny[k-2]*vec_nx[k-1])*acos(vec_nx[k-2]*vec_nx[k-1]+vec_ny[k-2]*vec_ny[k-1]);
							bnd_i.W2+=.5*int_k(0,0,vec_dt[k-1]);
							W_rs_nu dW2_10=get_W2_rs_type_II_k(1,0,{vec_x[k-2],vec_y[k-2]},{vec_x[k-1],vec_y[k-1]},{vec_nx[k-2],vec_ny[k-2]},vec_L[k-2],vec_dt[k-1]);
							bnd_i.W2_10[0]+=dW2_10.comps[0];
							bnd_i.W2_10[1]+=dW2_10.comps[1];
							W_rs_nu dW2_01=get_W2_rs_type_II_k(0,1,{vec_x[k-2],vec_y[k-2]},{vec_x[k-1],vec_y[k-1]},{vec_nx[k-2],vec_ny[k-2]},vec_L[k-2],vec_dt[k-1]);
							bnd_i.W2_01[0]+=dW2_01.comps[0];
							bnd_i.W2_01[1]+=dW2_01.comps[1];
						}
					}
					bnd_i.xi.push_back(vec_xi);
					bnd_i.dxi.push_back(vec_dxi);
					bnd_i.x.push_back(vec_x);
					bnd_i.y.push_back(vec_y);
					//
					bnd_i.L.push_back(vec_L);
					bnd_i.nx.push_back(vec_nx);
					bnd_i.ny.push_back(vec_ny);
					bnd_i.dt.push_back(vec_dt);
				}						
//				cout << "test 2b" << endl;
			}
			//
			// third case
			else if ((cmn_types[re_j(j,s)-1][1]==0)&(cmn_types[j][0]==1)) {
//				cout << "test 3a" << endl;
				int loc_cnt=1;
				if (thsr_ij[lorder[j]][0].size()==2) {
					if ((cmn_thsr[re_j(j,s)-1]<thsr_ij[lorder[j]][0][1][1])&(cmn_thsr[re_j(j,s)-1]>thsr_ij[lorder[j]][0][1][0])) {
						loc_cnt=2;
						bnd_i.th.push_back(linspace(cmn_thsr[re_j(j,s)-1],thsr_ij[lorder[j]][0][1][1],npts));
						bnd_i.th.push_back(linspace(thsr_ij[lorder[j]][0][0][0],thsr_ij[lorder[j]][0][0][1],npts));
					}
					else if ((cmn_thsr[re_j(j,s)-1]<thsr_ij[lorder[j]][0][0][1])&(cmn_thsr[re_j(j,s)-1]>thsr_ij[lorder[j]][0][0][0])) {
						bnd_i.th.push_back(linspace(cmn_thsr[re_j(j,s)-1],thsr_ij[lorder[j]][0][0][1],npts));
					}
				}
				else {
					if ((cmn_thsr[re_j(j,s)-1]<thsr_ij[lorder[j]][0][0][1])&(cmn_thsr[re_j(j,s)-1]>thsr_ij[lorder[j]][0][0][0])) {
						bnd_i.th.push_back(linspace(cmn_thsr[re_j(j,s)-1],thsr_ij[lorder[j]][0][0][1],npts));
					}
				}
				for (int i_cnt=loc_cnt;i_cnt>0;--i_cnt) {
					int ss=bnd_i.th.size();
					int n=bnd_i.th[ss-i_cnt].size();
					vector<double> vec_xi(n),vec_dxi(n),vec_x(n),vec_y(n);
					vector<double> vec_nx(n-1),vec_ny(n-1),vec_L(n-1),vec_dt(n-1);
					for (int k=0;k<n;++k) {
						vec_xi[k]=prm.xi_t0(bnd_i.th[ss-i_cnt][k]);
						vec_dxi[k]=prm.dxi_t0(bnd_i.th[ss-i_cnt][k]);
						vec_x[k]=prm.x_t0(bnd_i.th[ss-i_cnt][k]);
						vec_y[k]=prm.y_t0(bnd_i.th[ss-i_cnt][k]);
						//
						if (k>0) {
							vec_L[k-1]=pow(pow(vec_x[k]-vec_x[k-1],2.)+pow(vec_y[k]-vec_y[k-1],2.),.5);
							bnd_i.W1+=vec_L[k-1]/2.;
							vec_nx[k-1]=-(vec_y[k-1]-vec_y[k])/vec_L[k-1];
							vec_ny[k-1]=(vec_x[k-1]-vec_x[k])/vec_L[k-1];
							W_rs_nu dW1_10=get_W1_rs_type_II_k(1,0,{vec_x[k-1],vec_y[k-1]},{vec_x[k],vec_y[k]},{vec_nx[k-1],vec_ny[k-1]},vec_L[k-1]);
							bnd_i.W1_10[0]+=dW1_10.comps[0];
							bnd_i.W1_10[1]+=dW1_10.comps[1];	
							W_rs_nu dW1_01=get_W1_rs_type_II_k(0,1,{vec_x[k-1],vec_y[k-1]},{vec_x[k],vec_y[k]},{vec_nx[k-1],vec_ny[k-1]},vec_L[k-1]);
							bnd_i.W1_01[0]+=dW1_01.comps[0];
							bnd_i.W1_01[1]+=dW1_01.comps[1];
							bnd_i.W0+=1./4.*vec_L[k-1]*((vec_x[k-1]+vec_x[k])*vec_nx[k-1]+(vec_y[k-1]+vec_y[k])*vec_ny[k-1]);	
							W_rs_nu dW0_10=get_W0_rs_type_II_k(1,0,{vec_x[k-1],vec_y[k-1]},{vec_x[k],vec_y[k]},{vec_nx[k-1],vec_ny[k-1]},vec_L[k-1]);
							bnd_i.W0_10[0]+=dW0_10.comps[0];
							bnd_i.W0_10[1]+=dW0_10.comps[1];				
						}
						if ((k==1)&(bnd_i.L.size()>0)) {
							vec_dt[0]=sign(bnd_i.nx.back().back()*vec_ny[0]-bnd_i.ny.back().back()*vec_nx[0])*acos(bnd_i.nx.back().back()*vec_nx[0]+bnd_i.ny.back().back()*vec_ny[0]);
							bnd_i.W2+=.5*int_k(0,0,vec_dt[0]);
							W_rs_nu dW2_10=get_W2_rs_type_II_k(1,0,{bnd_i.x.back().end()[-2],bnd_i.y.back().end()[-2]},{bnd_i.x.back().back(),bnd_i.y.back().back()},{bnd_i.nx.back().back(),bnd_i.ny.back().back()},bnd_i.L.back().back(),vec_dt[k-1]);
							bnd_i.W2_10[0]+=dW2_10.comps[0];
							bnd_i.W2_10[1]+=dW2_10.comps[1];	
							W_rs_nu dW2_01=get_W2_rs_type_II_k(0,1,{bnd_i.x.back().end()[-2],bnd_i.y.back().end()[-2]},{bnd_i.x.back().back(),bnd_i.y.back().back()},{bnd_i.nx.back().back(),bnd_i.ny.back().back()},bnd_i.L.back().back(),vec_dt[k-1]);
							bnd_i.W2_01[0]+=dW2_01.comps[0];
							bnd_i.W2_01[1]+=dW2_01.comps[1];		
						} 
						else if (k>1) {
							vec_dt[k-1]=sign(vec_nx[k-2]*vec_ny[k-1]-vec_ny[k-2]*vec_nx[k-1])*acos(vec_nx[k-2]*vec_nx[k-1]+vec_ny[k-2]*vec_ny[k-1]);
							bnd_i.W2+=.5*int_k(0,0,vec_dt[k-1]);
							W_rs_nu dW2_10=get_W2_rs_type_II_k(1,0,{vec_x[k-2],vec_y[k-2]},{vec_x[k-1],vec_y[k-1]},{vec_nx[k-2],vec_ny[k-2]},vec_L[k-2],vec_dt[k-1]);
							bnd_i.W2_10[0]+=dW2_10.comps[0];
							bnd_i.W2_10[1]+=dW2_10.comps[1];
							W_rs_nu dW2_01=get_W2_rs_type_II_k(0,1,{vec_x[k-2],vec_y[k-2]},{vec_x[k-1],vec_y[k-1]},{vec_nx[k-2],vec_ny[k-2]},vec_L[k-2],vec_dt[k-1]);
							bnd_i.W2_01[0]+=dW2_01.comps[0];
							bnd_i.W2_01[1]+=dW2_01.comps[1];
						}
					}
					bnd_i.xi.push_back(vec_xi);
					bnd_i.dxi.push_back(vec_dxi);
					bnd_i.x.push_back(vec_x);
					bnd_i.y.push_back(vec_y);
					//
					bnd_i.L.push_back(vec_L);
					bnd_i.nx.push_back(vec_nx);
					bnd_i.ny.push_back(vec_ny);
					bnd_i.dt.push_back(vec_dt);
				}	
				
//				cout << "test 3b" << endl;					
				//
				loc_cnt=1;
				if (thsr_ij[lorder[j]][1].size()==2) {
					if ((cmn_thsr[j]<thsr_ij[lorder[j]][1][1][1])&(cmn_thsr[j]>thsr_ij[lorder[j]][1][1][0])) {
						loc_cnt=2;
						bnd_i.th.push_back(linspace(thsr_ij[lorder[j]][1][0][1],thsr_ij[lorder[j]][1][0][0],npts));
						bnd_i.th.push_back(linspace(thsr_ij[lorder[j]][1][1][1],cmn_thsr[j],npts));
					}
					else if ((cmn_thsr[j]<thsr_ij[lorder[j]][1][0][1])&(cmn_thsr[j]>thsr_ij[lorder[j]][1][0][0])) {
						bnd_i.th.push_back(linspace(thsr_ij[lorder[j]][1][0][1],cmn_thsr[j],npts));
					}
				}
				else {
					if ((cmn_thsr[j]<thsr_ij[lorder[j]][1][0][1])&(cmn_thsr[j]>thsr_ij[lorder[j]][1][0][0])) {
						bnd_i.th.push_back(linspace(thsr_ij[lorder[j]][1][0][1],cmn_thsr[j],npts));
					}
				}
				for (int i_cnt=loc_cnt;i_cnt>0;--i_cnt) {
					int ss=bnd_i.th.size();
					int n=bnd_i.th[ss-i_cnt].size();
					vector<double> vec_xi(n),vec_dxi(n),vec_x(n),vec_y(n);
					vector<double> vec_nx(n-1),vec_ny(n-1),vec_L(n-1),vec_dt(n-1);
					for (int k=0;k<n;++k) {
						vec_xi[k]=prm.xi_t1(bnd_i.th[ss-i_cnt][k]);
						vec_dxi[k]=prm.dxi_t1(bnd_i.th[ss-i_cnt][k]);
						vec_x[k]=prm.x_t1(bnd_i.th[ss-i_cnt][k]);
						vec_y[k]=prm.y_t1(bnd_i.th[ss-i_cnt][k]);
						//
						if (k>0) {
							vec_L[k-1]=pow(pow(vec_x[k]-vec_x[k-1],2.)+pow(vec_y[k]-vec_y[k-1],2.),.5);
							bnd_i.W1+=vec_L[k-1]/2.;
							vec_nx[k-1]=-(vec_y[k-1]-vec_y[k])/vec_L[k-1];
							vec_ny[k-1]=(vec_x[k-1]-vec_x[k])/vec_L[k-1];
							W_rs_nu dW1_10=get_W1_rs_type_II_k(1,0,{vec_x[k-1],vec_y[k-1]},{vec_x[k],vec_y[k]},{vec_nx[k-1],vec_ny[k-1]},vec_L[k-1]);
							bnd_i.W1_10[0]+=dW1_10.comps[0];
							bnd_i.W1_10[1]+=dW1_10.comps[1];	
							W_rs_nu dW1_01=get_W1_rs_type_II_k(0,1,{vec_x[k-1],vec_y[k-1]},{vec_x[k],vec_y[k]},{vec_nx[k-1],vec_ny[k-1]},vec_L[k-1]);
							bnd_i.W1_01[0]+=dW1_01.comps[0];
							bnd_i.W1_01[1]+=dW1_01.comps[1];
							bnd_i.W0+=1./4.*vec_L[k-1]*((vec_x[k-1]+vec_x[k])*vec_nx[k-1]+(vec_y[k-1]+vec_y[k])*vec_ny[k-1]);	
							W_rs_nu dW0_10=get_W0_rs_type_II_k(1,0,{vec_x[k-1],vec_y[k-1]},{vec_x[k],vec_y[k]},{vec_nx[k-1],vec_ny[k-1]},vec_L[k-1]);
							bnd_i.W0_10[0]+=dW0_10.comps[0];
							bnd_i.W0_10[1]+=dW0_10.comps[1];
							
							double B11,B22,B12;
							B11=vec_x[k-1]*vec_x[k-1]+vec_x[k]*vec_x[k]+vec_x[k-1]*vec_x[k];
							B22=vec_y[k-1]*vec_y[k-1]+vec_y[k]*vec_y[k]+vec_y[k-1]*vec_y[k];
							B12=vec_y[k-1]*vec_x[k-1]+vec_y[k]*vec_x[k]+.5*(vec_x[k-1]*vec_y[k]+vec_y[k-1]*vec_x[k]);
							bnd_i_W0_10[0]+=vec_L[k-1]/9.*(B11*vec_nx[k-1]+B12*vec_ny[k-1]);
							bnd_i_W0_10[1]+=vec_L[k-1]/9.*(B12*vec_nx[k-1]+B22*vec_ny[k-1]);
						}
						if ((k==1)&(bnd_i.L.size()>0)) {
							vec_dt[0]=sign(bnd_i.nx.back().back()*vec_ny[0]-bnd_i.ny.back().back()*vec_nx[0])*acos(bnd_i.nx.back().back()*vec_nx[0]+bnd_i.ny.back().back()*vec_ny[0]);
							bnd_i.W2+=.5*int_k(0,0,vec_dt[0]);
							W_rs_nu dW2_10=get_W2_rs_type_II_k(1,0,{bnd_i.x.back().end()[-2],bnd_i.y.back().end()[-2]},{bnd_i.x.back().back(),bnd_i.y.back().back()},{bnd_i.nx.back().back(),bnd_i.ny.back().back()},bnd_i.L.back().back(),vec_dt[k-1]);
							bnd_i.W2_10[0]+=dW2_10.comps[0];
							bnd_i.W2_10[1]+=dW2_10.comps[1];	
							W_rs_nu dW2_01=get_W2_rs_type_II_k(0,1,{bnd_i.x.back().end()[-2],bnd_i.y.back().end()[-2]},{bnd_i.x.back().back(),bnd_i.y.back().back()},{bnd_i.nx.back().back(),bnd_i.ny.back().back()},bnd_i.L.back().back(),vec_dt[k-1]);
							bnd_i.W2_01[0]+=dW2_01.comps[0];
							bnd_i.W2_01[1]+=dW2_01.comps[1];		
						} 
						else if (k>1) {
							vec_dt[k-1]=sign(vec_nx[k-2]*vec_ny[k-1]-vec_ny[k-2]*vec_nx[k-1])*acos(vec_nx[k-2]*vec_nx[k-1]+vec_ny[k-2]*vec_ny[k-1]);
							bnd_i.W2+=.5*int_k(0,0,vec_dt[k-1]);
							W_rs_nu dW2_10=get_W2_rs_type_II_k(1,0,{vec_x[k-2],vec_y[k-2]},{vec_x[k-1],vec_y[k-1]},{vec_nx[k-2],vec_ny[k-2]},vec_L[k-2],vec_dt[k-1]);
							bnd_i.W2_10[0]+=dW2_10.comps[0];
							bnd_i.W2_10[1]+=dW2_10.comps[1];
							W_rs_nu dW2_01=get_W2_rs_type_II_k(0,1,{vec_x[k-2],vec_y[k-2]},{vec_x[k-1],vec_y[k-1]},{vec_nx[k-2],vec_ny[k-2]},vec_L[k-2],vec_dt[k-1]);
							bnd_i.W2_01[0]+=dW2_01.comps[0];
							bnd_i.W2_01[1]+=dW2_01.comps[1];
						}
					}
					bnd_i.xi.push_back(vec_xi);
					bnd_i.dxi.push_back(vec_dxi);
					bnd_i.x.push_back(vec_x);
					bnd_i.y.push_back(vec_y);
					//
					bnd_i.L.push_back(vec_L);
					bnd_i.nx.push_back(vec_nx);
					bnd_i.ny.push_back(vec_ny);
					bnd_i.dt.push_back(vec_dt);
				}
//				cout << "test 3c" << endl;									
			}
		}
		bnd_i.dt[0][0]=sign(bnd_i.nx.back().back()*bnd_i.ny[0][0]-bnd_i.ny.back().back()*bnd_i.nx[0][0])*acos(bnd_i.nx.back().back()*bnd_i.nx[0][0]+bnd_i.ny.back().back()*bnd_i.ny[0][0]);
		//bnd_i.dt[0][0]=acos(bnd_i.nx.back().back()*bnd_i.nx[0][0]+bnd_i.ny.back().back()*bnd_i.ny[0][0]);
		bnd_i.W2+=.5*int_k(0,0,bnd_i.dt[0][0]);
		W_rs_nu dW2_10=get_W2_rs_type_II_k(1,0,{bnd_i.x.back().end()[-2],bnd_i.y.back().end()[-2]},{bnd_i.x[0][0],bnd_i.y[0][0]},{bnd_i.nx.back().back(),bnd_i.ny.back().back()},bnd_i.L.back().back(),bnd_i.dt[0][0]);
		bnd_i.W2_10[0]+=dW2_10.comps[0];
		bnd_i.W2_10[1]+=dW2_10.comps[1];
		//W_rs_nu dW2_01=get_W2_rs_type_II_k(0,1,{bnd_i.x.back().end()[-2],bnd_i.y.back().end()[-2]},{bnd_i.x[0][0],bnd_i.y[0][0]},{bnd_i.nx.back().back(),bnd_i.ny.back().back()},bnd_i.L.back().back(),bnd_i.dt[0][0]);
		W_rs_nu dW2_01=get_W2_rs_type_II_k(0,1,{bnd_i.x.back().end()[-2],bnd_i.y.back().end()[-2]},{bnd_i.x.back().back(),bnd_i.y.back().back()},{bnd_i.nx.back().back(),bnd_i.ny.back().back()},bnd_i.L.back().back(),bnd_i.dt[0][0]);
		bnd_i.W2_01[0]+=dW2_01.comps[0];
		bnd_i.W2_01[1]+=dW2_01.comps[1];
		//
		bnd_i.W0_10[0]/=bnd_i.W0;
		bnd_i.W0_10[1]/=bnd_i.W0;
		bnd_i.W1_10[0]/=bnd_i.W1;
		bnd_i.W1_10[1]/=bnd_i.W1;
		bnd_i.W2_10[0]/=bnd_i.W2;
		bnd_i.W2_10[1]/=bnd_i.W2;
		
		for (size_t j=0;j<bnd_i.x.size();++j) {
			bnd_i.x_w0.push_back(bnd_i.x[j]);
			bnd_i.y_w0.push_back(bnd_i.y[j]);
			bnd_i.x_w1.push_back(bnd_i.x[j]);
			bnd_i.y_w1.push_back(bnd_i.y[j]);
			bnd_i.x_w2.push_back(bnd_i.x[j]);
			bnd_i.y_w2.push_back(bnd_i.y[j]);
			for (size_t k=0;k<bnd_i.x[j].size();++k) {
				bnd_i.x_w0[j][k]-=bnd_i.W0_10[0];
				bnd_i.y_w0[j][k]-=bnd_i.W0_10[1];
				bnd_i.x_w1[j][k]-=bnd_i.W1_10[0];
				bnd_i.y_w1[j][k]-=bnd_i.W1_10[1];
				bnd_i.x_w2[j][k]-=bnd_i.W2_10[0];
				bnd_i.y_w2[j][k]-=bnd_i.W2_10[1];
			}
		}		
		//
//		cout << "test a" << endl;
		for (size_t j=0;j<bnd_i.x.size();++j) {
			for (size_t k=0;k<bnd_i.x[j].size();++k) {
				if (k>0) {		
					//
					// Tensors of type 0, w0_rs with r+s<=6
					// order 2
					W_rs_nu dW0_20=get_W0_rs_type_II_k(2,0,{bnd_i.x_w0[j][k-1],bnd_i.y_w0[j][k-1]},{bnd_i.x_w0[j][k],bnd_i.y_w0[j][k]},{bnd_i.nx[j][k-1],bnd_i.ny[j][k-1]},bnd_i.L[j][k-1]);
					bnd_i.W0_20[0]+=dW0_20.comps[0];
					bnd_i.W0_20[1]+=dW0_20.comps[1];
					bnd_i.W0_20[2]+=dW0_20.comps[2];
					// order 3
					W_rs_nu dW0_30=get_W0_rs_type_II_k(3,0,{bnd_i.x_w0[j][k-1],bnd_i.y_w0[j][k-1]},{bnd_i.x_w0[j][k],bnd_i.y_w0[j][k]},{bnd_i.nx[j][k-1],bnd_i.ny[j][k-1]},bnd_i.L[j][k-1]);
					bnd_i.W0_30[0]+=dW0_30.comps[0];
					bnd_i.W0_30[1]+=dW0_30.comps[1];
					bnd_i.W0_30[2]+=dW0_30.comps[2];
					bnd_i.W0_30[3]+=dW0_30.comps[3];
					// order 4
					W_rs_nu dW0_40=get_W0_rs_type_II_k(4,0,{bnd_i.x_w0[j][k-1],bnd_i.y_w0[j][k-1]},{bnd_i.x_w0[j][k],bnd_i.y_w0[j][k]},{bnd_i.nx[j][k-1],bnd_i.ny[j][k-1]},bnd_i.L[j][k-1]);
					bnd_i.W0_40[0]+=dW0_40.comps[0];
					bnd_i.W0_40[1]+=dW0_40.comps[1];
					bnd_i.W0_40[2]+=dW0_40.comps[2];
					bnd_i.W0_40[3]+=dW0_40.comps[3];
					bnd_i.W0_40[4]+=dW0_40.comps[4];
					// order 5
					W_rs_nu dW0_50=get_W0_rs_type_II_k(5,0,{bnd_i.x_w0[j][k-1],bnd_i.y_w0[j][k-1]},{bnd_i.x_w0[j][k],bnd_i.y_w0[j][k]},{bnd_i.nx[j][k-1],bnd_i.ny[j][k-1]},bnd_i.L[j][k-1]);
					bnd_i.W0_50[0]+=dW0_50.comps[0];
					bnd_i.W0_50[1]+=dW0_50.comps[1];
					bnd_i.W0_50[2]+=dW0_50.comps[2];
					bnd_i.W0_50[3]+=dW0_50.comps[3];
					bnd_i.W0_50[4]+=dW0_50.comps[4];
					bnd_i.W0_50[5]+=dW0_50.comps[5];
					// order 6
					W_rs_nu dW0_60=get_W0_rs_type_II_k(6,0,{bnd_i.x_w0[j][k-1],bnd_i.y_w0[j][k-1]},{bnd_i.x_w0[j][k],bnd_i.y_w0[j][k]},{bnd_i.nx[j][k-1],bnd_i.ny[j][k-1]},bnd_i.L[j][k-1]);
					bnd_i.W0_60[0]+=dW0_60.comps[0];
					bnd_i.W0_60[1]+=dW0_60.comps[1];
					bnd_i.W0_60[2]+=dW0_60.comps[2];
					bnd_i.W0_60[3]+=dW0_60.comps[3];
					bnd_i.W0_60[4]+=dW0_60.comps[4];
					bnd_i.W0_60[5]+=dW0_60.comps[5];
					bnd_i.W0_60[6]+=dW0_60.comps[6];
					//
					// Tensors of type I, w1_rs with r+s<=6
					// order 2
					W_rs_nu dW1_20=get_W1_rs_type_II_k(2,0,{bnd_i.x_w1[j][k-1],bnd_i.y_w1[j][k-1]},{bnd_i.x_w1[j][k],bnd_i.y_w1[j][k]},{bnd_i.nx[j][k-1],bnd_i.ny[j][k-1]},bnd_i.L[j][k-1]);
					bnd_i.W1_20[0]+=dW1_20.comps[0];
					bnd_i.W1_20[1]+=dW1_20.comps[1];
					bnd_i.W1_20[2]+=dW1_20.comps[2];
					W_rs_nu dW1_11=get_W1_rs_type_II_k(1,1,{bnd_i.x_w1[j][k-1],bnd_i.y_w1[j][k-1]},{bnd_i.x_w1[j][k],bnd_i.y_w1[j][k]},{bnd_i.nx[j][k-1],bnd_i.ny[j][k-1]},bnd_i.L[j][k-1]);
					bnd_i.W1_11[0]+=dW1_11.comps[0];
					bnd_i.W1_11[1]+=dW1_11.comps[1];
					bnd_i.W1_11[2]+=dW1_11.comps[2];
					W_rs_nu dW1_02=get_W1_rs_type_II_k(0,2,{bnd_i.x_w1[j][k-1],bnd_i.y_w1[j][k-1]},{bnd_i.x_w1[j][k],bnd_i.y_w1[j][k]},{bnd_i.nx[j][k-1],bnd_i.ny[j][k-1]},bnd_i.L[j][k-1]);
					bnd_i.W1_02[0]+=dW1_02.comps[0];
					bnd_i.W1_02[1]+=dW1_02.comps[1];
					bnd_i.W1_02[2]+=dW1_02.comps[2];
					// order 3
					W_rs_nu dW1_30=get_W1_rs_type_II_k(3,0,{bnd_i.x_w1[j][k-1],bnd_i.y_w1[j][k-1]},{bnd_i.x_w1[j][k],bnd_i.y_w1[j][k]},{bnd_i.nx[j][k-1],bnd_i.ny[j][k-1]},bnd_i.L[j][k-1]);
					bnd_i.W1_30[0]+=dW1_30.comps[0];
					bnd_i.W1_30[1]+=dW1_30.comps[1];
					bnd_i.W1_30[2]+=dW1_30.comps[2];
					bnd_i.W1_30[3]+=dW1_30.comps[3];
					W_rs_nu dW1_21=get_W1_rs_type_II_k(2,1,{bnd_i.x_w1[j][k-1],bnd_i.y_w1[j][k-1]},{bnd_i.x_w1[j][k],bnd_i.y_w1[j][k]},{bnd_i.nx[j][k-1],bnd_i.ny[j][k-1]},bnd_i.L[j][k-1]);
					bnd_i.W1_21[0]+=dW1_21.comps[0];
					bnd_i.W1_21[1]+=dW1_21.comps[1];
					bnd_i.W1_21[2]+=dW1_21.comps[2];
					bnd_i.W1_21[3]+=dW1_21.comps[3];
					W_rs_nu dW1_12=get_W1_rs_type_II_k(1,2,{bnd_i.x_w1[j][k-1],bnd_i.y_w1[j][k-1]},{bnd_i.x_w1[j][k],bnd_i.y_w1[j][k]},{bnd_i.nx[j][k-1],bnd_i.ny[j][k-1]},bnd_i.L[j][k-1]);
					bnd_i.W1_12[0]+=dW1_12.comps[0];
					bnd_i.W1_12[1]+=dW1_12.comps[1];
					bnd_i.W1_12[2]+=dW1_12.comps[2];
					bnd_i.W1_12[3]+=dW1_12.comps[3];
					W_rs_nu dW1_03=get_W1_rs_type_II_k(0,3,{bnd_i.x_w1[j][k-1],bnd_i.y_w1[j][k-1]},{bnd_i.x_w1[j][k],bnd_i.y_w1[j][k]},{bnd_i.nx[j][k-1],bnd_i.ny[j][k-1]},bnd_i.L[j][k-1]);
					bnd_i.W1_03[0]+=dW1_03.comps[0];
					bnd_i.W1_03[1]+=dW1_03.comps[1];
					bnd_i.W1_03[2]+=dW1_03.comps[2];
					bnd_i.W1_03[3]+=dW1_03.comps[3];
					// order 4
					W_rs_nu dW1_40=get_W1_rs_type_II_k(4,0,{bnd_i.x_w1[j][k-1],bnd_i.y_w1[j][k-1]},{bnd_i.x_w1[j][k],bnd_i.y_w1[j][k]},{bnd_i.nx[j][k-1],bnd_i.ny[j][k-1]},bnd_i.L[j][k-1]);
					bnd_i.W1_40[0]+=dW1_40.comps[0];
					bnd_i.W1_40[1]+=dW1_40.comps[1];
					bnd_i.W1_40[2]+=dW1_40.comps[2];
					bnd_i.W1_40[3]+=dW1_40.comps[3];
					bnd_i.W1_40[4]+=dW1_40.comps[4];
					W_rs_nu dW1_31=get_W1_rs_type_II_k(3,1,{bnd_i.x_w1[j][k-1],bnd_i.y_w1[j][k-1]},{bnd_i.x_w1[j][k],bnd_i.y_w1[j][k]},{bnd_i.nx[j][k-1],bnd_i.ny[j][k-1]},bnd_i.L[j][k-1]);
					bnd_i.W1_31[0]+=dW1_31.comps[0];
					bnd_i.W1_31[1]+=dW1_31.comps[1];
					bnd_i.W1_31[2]+=dW1_31.comps[2];
					bnd_i.W1_31[3]+=dW1_31.comps[3];
					bnd_i.W1_31[4]+=dW1_31.comps[4];
					W_rs_nu dW1_22=get_W1_rs_type_II_k(2,2,{bnd_i.x_w1[j][k-1],bnd_i.y_w1[j][k-1]},{bnd_i.x_w1[j][k],bnd_i.y_w1[j][k]},{bnd_i.nx[j][k-1],bnd_i.ny[j][k-1]},bnd_i.L[j][k-1]);
					bnd_i.W1_22[0]+=dW1_22.comps[0];
					bnd_i.W1_22[1]+=dW1_22.comps[1];
					bnd_i.W1_22[2]+=dW1_22.comps[2];
					bnd_i.W1_22[3]+=dW1_22.comps[3];
					bnd_i.W1_22[4]+=dW1_22.comps[4];
					W_rs_nu dW1_13=get_W1_rs_type_II_k(1,3,{bnd_i.x_w1[j][k-1],bnd_i.y_w1[j][k-1]},{bnd_i.x_w1[j][k],bnd_i.y_w1[j][k]},{bnd_i.nx[j][k-1],bnd_i.ny[j][k-1]},bnd_i.L[j][k-1]);
					bnd_i.W1_13[0]+=dW1_13.comps[0];
					bnd_i.W1_13[1]+=dW1_13.comps[1];
					bnd_i.W1_13[2]+=dW1_13.comps[2];
					bnd_i.W1_13[3]+=dW1_13.comps[3];
					bnd_i.W1_13[4]+=dW1_13.comps[4];
					W_rs_nu dW1_04=get_W1_rs_type_II_k(0,4,{bnd_i.x_w1[j][k-1],bnd_i.y_w1[j][k-1]},{bnd_i.x_w1[j][k],bnd_i.y_w1[j][k]},{bnd_i.nx[j][k-1],bnd_i.ny[j][k-1]},bnd_i.L[j][k-1]);
					bnd_i.W1_04[0]+=dW1_04.comps[0];
					bnd_i.W1_04[1]+=dW1_04.comps[1];
					bnd_i.W1_04[2]+=dW1_04.comps[2];
					bnd_i.W1_04[3]+=dW1_04.comps[3];
					bnd_i.W1_04[4]+=dW1_04.comps[4];
					// order 5
					W_rs_nu dW1_50=get_W1_rs_type_II_k(5,0,{bnd_i.x_w1[j][k-1],bnd_i.y_w1[j][k-1]},{bnd_i.x_w1[j][k],bnd_i.y_w1[j][k]},{bnd_i.nx[j][k-1],bnd_i.ny[j][k-1]},bnd_i.L[j][k-1]);
					bnd_i.W1_50[0]+=dW1_50.comps[0];
					bnd_i.W1_50[1]+=dW1_50.comps[1];
					bnd_i.W1_50[2]+=dW1_50.comps[2];
					bnd_i.W1_50[3]+=dW1_50.comps[3];
					bnd_i.W1_50[4]+=dW1_50.comps[4];
					bnd_i.W1_50[5]+=dW1_50.comps[5];
					W_rs_nu dW1_41=get_W1_rs_type_II_k(4,1,{bnd_i.x_w1[j][k-1],bnd_i.y_w1[j][k-1]},{bnd_i.x_w1[j][k],bnd_i.y_w1[j][k]},{bnd_i.nx[j][k-1],bnd_i.ny[j][k-1]},bnd_i.L[j][k-1]);
					bnd_i.W1_41[0]+=dW1_41.comps[0];
					bnd_i.W1_41[1]+=dW1_41.comps[1];
					bnd_i.W1_41[2]+=dW1_41.comps[2];
					bnd_i.W1_41[3]+=dW1_41.comps[3];
					bnd_i.W1_41[4]+=dW1_41.comps[4];
					bnd_i.W1_41[5]+=dW1_41.comps[5];
					W_rs_nu dW1_32=get_W1_rs_type_II_k(3,2,{bnd_i.x_w1[j][k-1],bnd_i.y_w1[j][k-1]},{bnd_i.x_w1[j][k],bnd_i.y_w1[j][k]},{bnd_i.nx[j][k-1],bnd_i.ny[j][k-1]},bnd_i.L[j][k-1]);
					bnd_i.W1_32[0]+=dW1_32.comps[0];
					bnd_i.W1_32[1]+=dW1_32.comps[1];
					bnd_i.W1_32[2]+=dW1_32.comps[2];
					bnd_i.W1_32[3]+=dW1_32.comps[3];
					bnd_i.W1_32[4]+=dW1_32.comps[4];
					bnd_i.W1_32[5]+=dW1_32.comps[5];
					W_rs_nu dW1_23=get_W1_rs_type_II_k(2,3,{bnd_i.x_w1[j][k-1],bnd_i.y_w1[j][k-1]},{bnd_i.x_w1[j][k],bnd_i.y_w1[j][k]},{bnd_i.nx[j][k-1],bnd_i.ny[j][k-1]},bnd_i.L[j][k-1]);
					bnd_i.W1_23[0]+=dW1_23.comps[0];
					bnd_i.W1_23[1]+=dW1_23.comps[1];
					bnd_i.W1_23[2]+=dW1_23.comps[2];
					bnd_i.W1_23[3]+=dW1_23.comps[3];
					bnd_i.W1_23[4]+=dW1_23.comps[4];
					bnd_i.W1_23[5]+=dW1_23.comps[5];
					W_rs_nu dW1_14=get_W1_rs_type_II_k(1,4,{bnd_i.x_w1[j][k-1],bnd_i.y_w1[j][k-1]},{bnd_i.x_w1[j][k],bnd_i.y_w1[j][k]},{bnd_i.nx[j][k-1],bnd_i.ny[j][k-1]},bnd_i.L[j][k-1]);
					bnd_i.W1_14[0]+=dW1_14.comps[0];
					bnd_i.W1_14[1]+=dW1_14.comps[1];
					bnd_i.W1_14[2]+=dW1_14.comps[2];
					bnd_i.W1_14[3]+=dW1_14.comps[3];
					bnd_i.W1_14[4]+=dW1_14.comps[4];
					bnd_i.W1_14[5]+=dW1_14.comps[5];
					W_rs_nu dW1_05=get_W1_rs_type_II_k(0,5,{bnd_i.x_w1[j][k-1],bnd_i.y_w1[j][k-1]},{bnd_i.x_w1[j][k],bnd_i.y_w1[j][k]},{bnd_i.nx[j][k-1],bnd_i.ny[j][k-1]},bnd_i.L[j][k-1]);
					bnd_i.W1_05[0]+=dW1_05.comps[0];
					bnd_i.W1_05[1]+=dW1_05.comps[1];
					bnd_i.W1_05[2]+=dW1_05.comps[2];
					bnd_i.W1_05[3]+=dW1_05.comps[3];
					bnd_i.W1_05[4]+=dW1_05.comps[4];
					bnd_i.W1_05[5]+=dW1_05.comps[5];
					// order 6
					W_rs_nu dW1_60=get_W1_rs_type_II_k(6,0,{bnd_i.x_w1[j][k-1],bnd_i.y_w1[j][k-1]},{bnd_i.x_w1[j][k],bnd_i.y_w1[j][k]},{bnd_i.nx[j][k-1],bnd_i.ny[j][k-1]},bnd_i.L[j][k-1]);
					bnd_i.W1_60[0]+=dW1_60.comps[0];
					bnd_i.W1_60[1]+=dW1_60.comps[1];
					bnd_i.W1_60[2]+=dW1_60.comps[2];
					bnd_i.W1_60[3]+=dW1_60.comps[3];
					bnd_i.W1_60[4]+=dW1_60.comps[4];
					bnd_i.W1_60[5]+=dW1_60.comps[5];
					bnd_i.W1_60[6]+=dW1_60.comps[6];
					W_rs_nu dW1_51=get_W1_rs_type_II_k(5,1,{bnd_i.x_w1[j][k-1],bnd_i.y_w1[j][k-1]},{bnd_i.x_w1[j][k],bnd_i.y_w1[j][k]},{bnd_i.nx[j][k-1],bnd_i.ny[j][k-1]},bnd_i.L[j][k-1]);
					bnd_i.W1_51[0]+=dW1_51.comps[0];
					bnd_i.W1_51[1]+=dW1_51.comps[1];
					bnd_i.W1_51[2]+=dW1_51.comps[2];
					bnd_i.W1_51[3]+=dW1_51.comps[3];
					bnd_i.W1_51[4]+=dW1_51.comps[4];
					bnd_i.W1_51[5]+=dW1_51.comps[5];
					bnd_i.W1_51[6]+=dW1_51.comps[6];
					W_rs_nu dW1_42=get_W1_rs_type_II_k(4,2,{bnd_i.x_w1[j][k-1],bnd_i.y_w1[j][k-1]},{bnd_i.x_w1[j][k],bnd_i.y_w1[j][k]},{bnd_i.nx[j][k-1],bnd_i.ny[j][k-1]},bnd_i.L[j][k-1]);
					bnd_i.W1_42[0]+=dW1_42.comps[0];
					bnd_i.W1_42[1]+=dW1_42.comps[1];
					bnd_i.W1_42[2]+=dW1_42.comps[2];
					bnd_i.W1_42[3]+=dW1_42.comps[3];
					bnd_i.W1_42[4]+=dW1_42.comps[4];
					bnd_i.W1_42[5]+=dW1_42.comps[5];
					bnd_i.W1_42[6]+=dW1_42.comps[6];
					W_rs_nu dW1_33=get_W1_rs_type_II_k(3,3,{bnd_i.x_w1[j][k-1],bnd_i.y_w1[j][k-1]},{bnd_i.x_w1[j][k],bnd_i.y_w1[j][k]},{bnd_i.nx[j][k-1],bnd_i.ny[j][k-1]},bnd_i.L[j][k-1]);
					bnd_i.W1_33[0]+=dW1_33.comps[0];
					bnd_i.W1_33[1]+=dW1_33.comps[1];
					bnd_i.W1_33[2]+=dW1_33.comps[2];
					bnd_i.W1_33[3]+=dW1_33.comps[3];
					bnd_i.W1_33[4]+=dW1_33.comps[4];
					bnd_i.W1_33[5]+=dW1_33.comps[5];
					bnd_i.W1_33[6]+=dW1_33.comps[6];
					W_rs_nu dW1_24=get_W1_rs_type_II_k(2,4,{bnd_i.x_w1[j][k-1],bnd_i.y_w1[j][k-1]},{bnd_i.x_w1[j][k],bnd_i.y_w1[j][k]},{bnd_i.nx[j][k-1],bnd_i.ny[j][k-1]},bnd_i.L[j][k-1]);
					bnd_i.W1_24[0]+=dW1_24.comps[0];
					bnd_i.W1_24[1]+=dW1_24.comps[1];
					bnd_i.W1_24[2]+=dW1_24.comps[2];
					bnd_i.W1_24[3]+=dW1_24.comps[3];
					bnd_i.W1_24[4]+=dW1_24.comps[4];
					bnd_i.W1_24[5]+=dW1_24.comps[5];
					bnd_i.W1_24[6]+=dW1_24.comps[6];
					W_rs_nu dW1_15=get_W1_rs_type_II_k(1,5,{bnd_i.x_w1[j][k-1],bnd_i.y_w1[j][k-1]},{bnd_i.x_w1[j][k],bnd_i.y_w1[j][k]},{bnd_i.nx[j][k-1],bnd_i.ny[j][k-1]},bnd_i.L[j][k-1]);
					bnd_i.W1_15[0]+=dW1_15.comps[0];
					bnd_i.W1_15[1]+=dW1_15.comps[1];
					bnd_i.W1_15[2]+=dW1_15.comps[2];
					bnd_i.W1_15[3]+=dW1_15.comps[3];
					bnd_i.W1_15[4]+=dW1_15.comps[4];
					bnd_i.W1_15[5]+=dW1_15.comps[5];
					bnd_i.W1_15[6]+=dW1_15.comps[6];
					W_rs_nu dW1_06=get_W1_rs_type_II_k(0,6,{bnd_i.x_w1[j][k-1],bnd_i.y_w1[j][k-1]},{bnd_i.x_w1[j][k],bnd_i.y_w1[j][k]},{bnd_i.nx[j][k-1],bnd_i.ny[j][k-1]},bnd_i.L[j][k-1]);
					bnd_i.W1_06[0]+=dW1_06.comps[0];
					bnd_i.W1_06[1]+=dW1_06.comps[1];
					bnd_i.W1_06[2]+=dW1_06.comps[2];
					bnd_i.W1_06[3]+=dW1_06.comps[3];
					bnd_i.W1_06[4]+=dW1_06.comps[4];
					bnd_i.W1_06[5]+=dW1_06.comps[5];
					bnd_i.W1_06[6]+=dW1_06.comps[6];
				}
				//
				// Tensors of type II, w2_rs with r+s<=6
				// order 2
				vector<double> nk(2),xk(2),xk_p1(2);
				double Lk=0.;
				//
				if ((k==0)&(j>0)) {
					nk={bnd_i.nx[j-1].back(),bnd_i.ny[j-1].back()};
					Lk=bnd_i.L[j-1].back();
					xk={bnd_i.x_w2[j-1].end()[-2],bnd_i.y_w2[j-1].end()[-2]};
					xk_p1={bnd_i.x_w2[j-1].back(),bnd_i.y_w2[j-1].back()};
				}
				else if ((k==0)&(j==0)) {
					nk={bnd_i.nx.back().back(),bnd_i.ny.back().back()};
					Lk=bnd_i.L.back().back();
					xk={bnd_i.x_w2.back().end()[-2],bnd_i.y_w2.back().end()[-2]};
					xk_p1={bnd_i.x_w2.back().back(),bnd_i.y_w2.back().back()};
				}
				else {					
					nk={bnd_i.nx[j][k-1],bnd_i.ny[j][k-1]};					
					Lk=bnd_i.L[j][k-1];
					xk={bnd_i.x_w2[j][k-1],bnd_i.y_w2[j][k-1]};
					xk_p1={bnd_i.x_w2[j][k],bnd_i.y_w2[j][k]};
				}
				W_rs_nu dW2_20=get_W2_rs_type_II_k(2,0,xk,xk_p1,nk,Lk,bnd_i.dt[j][k]);	
				bnd_i.W2_20[0]+=dW2_20.comps[0];
				bnd_i.W2_20[1]+=dW2_20.comps[1];
				bnd_i.W2_20[2]+=dW2_20.comps[2];
				W_rs_nu dW2_11=get_W2_rs_type_II_k(1,1,xk,xk_p1,nk,Lk,bnd_i.dt[j][k]);	
				bnd_i.W2_11[0]+=dW2_11.comps[0];
				bnd_i.W2_11[1]+=dW2_11.comps[1];
				bnd_i.W2_11[2]+=dW2_11.comps[2];
				W_rs_nu dW2_02=get_W2_rs_type_II_k(0,2,xk,xk_p1,nk,Lk,bnd_i.dt[j][k]);
				bnd_i.W2_02[0]+=dW2_02.comps[0];
				bnd_i.W2_02[1]+=dW2_02.comps[1];
				bnd_i.W2_02[2]+=dW2_02.comps[2];
				// order 3
				W_rs_nu dW2_30=get_W2_rs_type_II_k(3,0,xk,xk_p1,nk,Lk,bnd_i.dt[j][k]);
				bnd_i.W2_30[0]+=dW2_30.comps[0];
				bnd_i.W2_30[1]+=dW2_30.comps[1];
				bnd_i.W2_30[2]+=dW2_30.comps[2];
				bnd_i.W2_30[3]+=dW2_30.comps[3];
				W_rs_nu dW2_21=get_W2_rs_type_II_k(2,1,xk,xk_p1,nk,Lk,bnd_i.dt[j][k]);
				bnd_i.W2_21[0]+=dW2_21.comps[0];
				bnd_i.W2_21[1]+=dW2_21.comps[1];
				bnd_i.W2_21[2]+=dW2_21.comps[2];
				bnd_i.W2_21[3]+=dW2_21.comps[3];
				W_rs_nu dW2_12=get_W2_rs_type_II_k(1,2,xk,xk_p1,nk,Lk,bnd_i.dt[j][k]);
				bnd_i.W2_12[0]+=dW2_12.comps[0];
				bnd_i.W2_12[1]+=dW2_12.comps[1];
				bnd_i.W2_12[2]+=dW2_12.comps[2];
				bnd_i.W2_12[3]+=dW2_12.comps[3];
				W_rs_nu dW2_03=get_W2_rs_type_II_k(0,3,xk,xk_p1,nk,Lk,bnd_i.dt[j][k]);
				bnd_i.W2_03[0]+=dW2_03.comps[0];
				bnd_i.W2_03[1]+=dW2_03.comps[1];
				bnd_i.W2_03[2]+=dW2_03.comps[2];
				bnd_i.W2_03[3]+=dW2_03.comps[3];
				// order 4
				W_rs_nu dW2_40=get_W2_rs_type_II_k(4,0,xk,xk_p1,nk,Lk,bnd_i.dt[j][k]);
				bnd_i.W2_40[0]+=dW2_40.comps[0];
				bnd_i.W2_40[1]+=dW2_40.comps[1];
				bnd_i.W2_40[2]+=dW2_40.comps[2];
				bnd_i.W2_40[3]+=dW2_40.comps[3];
				bnd_i.W2_40[4]+=dW2_40.comps[4];
				W_rs_nu dW2_31=get_W2_rs_type_II_k(3,1,xk,xk_p1,nk,Lk,bnd_i.dt[j][k]);
				bnd_i.W2_31[0]+=dW2_31.comps[0];
				bnd_i.W2_31[1]+=dW2_31.comps[1];
				bnd_i.W2_31[2]+=dW2_31.comps[2];
				bnd_i.W2_31[3]+=dW2_31.comps[3];
				bnd_i.W2_31[4]+=dW2_31.comps[4];
				W_rs_nu dW2_22=get_W2_rs_type_II_k(2,2,xk,xk_p1,nk,Lk,bnd_i.dt[j][k]);
				bnd_i.W2_22[0]+=dW2_22.comps[0];
				bnd_i.W2_22[1]+=dW2_22.comps[1];
				bnd_i.W2_22[2]+=dW2_22.comps[2];
				bnd_i.W2_22[3]+=dW2_22.comps[3];
				bnd_i.W2_22[4]+=dW2_22.comps[4];
				W_rs_nu dW2_13=get_W2_rs_type_II_k(1,3,xk,xk_p1,nk,Lk,bnd_i.dt[j][k]);
				bnd_i.W2_13[0]+=dW2_13.comps[0];
				bnd_i.W2_13[1]+=dW2_13.comps[1];
				bnd_i.W2_13[2]+=dW2_13.comps[2];
				bnd_i.W2_13[3]+=dW2_13.comps[3];
				bnd_i.W2_13[4]+=dW2_13.comps[4];
				W_rs_nu dW2_04=get_W2_rs_type_II_k(0,4,xk,xk_p1,nk,Lk,bnd_i.dt[j][k]);
				bnd_i.W2_04[0]+=dW2_04.comps[0];
				bnd_i.W2_04[1]+=dW2_04.comps[1];
				bnd_i.W2_04[2]+=dW2_04.comps[2];
				bnd_i.W2_04[3]+=dW2_04.comps[3];
				bnd_i.W2_04[4]+=dW2_04.comps[4];
				// order 5
				W_rs_nu dW2_50=get_W2_rs_type_II_k(5,0,xk,xk_p1,nk,Lk,bnd_i.dt[j][k]);
				bnd_i.W2_50[0]+=dW2_50.comps[0];
				bnd_i.W2_50[1]+=dW2_50.comps[1];
				bnd_i.W2_50[2]+=dW2_50.comps[2];
				bnd_i.W2_50[3]+=dW2_50.comps[3];
				bnd_i.W2_50[4]+=dW2_50.comps[4];
				bnd_i.W2_50[5]+=dW2_50.comps[5];
				W_rs_nu dW2_41=get_W2_rs_type_II_k(4,1,xk,xk_p1,nk,Lk,bnd_i.dt[j][k]);
				bnd_i.W2_41[0]+=dW2_41.comps[0];
				bnd_i.W2_41[1]+=dW2_41.comps[1];
				bnd_i.W2_41[2]+=dW2_41.comps[2];
				bnd_i.W2_41[3]+=dW2_41.comps[3];
				bnd_i.W2_41[4]+=dW2_41.comps[4];
				bnd_i.W2_41[5]+=dW2_41.comps[5];
				W_rs_nu dW2_32=get_W2_rs_type_II_k(3,2,xk,xk_p1,nk,Lk,bnd_i.dt[j][k]);
				bnd_i.W2_32[0]+=dW2_32.comps[0];
				bnd_i.W2_32[1]+=dW2_32.comps[1];
				bnd_i.W2_32[2]+=dW2_32.comps[2];
				bnd_i.W2_32[3]+=dW2_32.comps[3];
				bnd_i.W2_32[4]+=dW2_32.comps[4];
				bnd_i.W2_32[5]+=dW2_32.comps[5];
				W_rs_nu dW2_23=get_W2_rs_type_II_k(2,3,xk,xk_p1,nk,Lk,bnd_i.dt[j][k]);
				bnd_i.W2_23[0]+=dW2_23.comps[0];
				bnd_i.W2_23[1]+=dW2_23.comps[1];
				bnd_i.W2_23[2]+=dW2_23.comps[2];
				bnd_i.W2_23[3]+=dW2_23.comps[3];
				bnd_i.W2_23[4]+=dW2_23.comps[4];
				bnd_i.W2_23[5]+=dW2_23.comps[5];
				W_rs_nu dW2_14=get_W2_rs_type_II_k(1,4,xk,xk_p1,nk,Lk,bnd_i.dt[j][k]);
				bnd_i.W2_14[0]+=dW2_14.comps[0];
				bnd_i.W2_14[1]+=dW2_14.comps[1];
				bnd_i.W2_14[2]+=dW2_14.comps[2];
				bnd_i.W2_14[3]+=dW2_14.comps[3];
				bnd_i.W2_14[4]+=dW2_14.comps[4];
				bnd_i.W2_14[5]+=dW2_14.comps[5];
				W_rs_nu dW2_05=get_W2_rs_type_II_k(0,5,xk,xk_p1,nk,Lk,bnd_i.dt[j][k]);
				bnd_i.W2_05[0]+=dW2_05.comps[0];
				bnd_i.W2_05[1]+=dW2_05.comps[1];
				bnd_i.W2_05[2]+=dW2_05.comps[2];
				bnd_i.W2_05[3]+=dW2_05.comps[3];
				bnd_i.W2_05[4]+=dW2_05.comps[4];
				bnd_i.W2_05[5]+=dW2_05.comps[5];
				// order 6
				W_rs_nu dW2_60=get_W2_rs_type_II_k(6,0,xk,xk_p1,nk,Lk,bnd_i.dt[j][k]);
				bnd_i.W2_60[0]+=dW2_60.comps[0];
				bnd_i.W2_60[1]+=dW2_60.comps[1];
				bnd_i.W2_60[2]+=dW2_60.comps[2];
				bnd_i.W2_60[3]+=dW2_60.comps[3];
				bnd_i.W2_60[4]+=dW2_60.comps[4];
				bnd_i.W2_60[5]+=dW2_60.comps[5];
				bnd_i.W2_60[6]+=dW2_60.comps[6];
				W_rs_nu dW2_51=get_W2_rs_type_II_k(5,1,xk,xk_p1,nk,Lk,bnd_i.dt[j][k]);
				bnd_i.W2_51[0]+=dW2_51.comps[0];
				bnd_i.W2_51[1]+=dW2_51.comps[1];
				bnd_i.W2_51[2]+=dW2_51.comps[2];
				bnd_i.W2_51[3]+=dW2_51.comps[3];
				bnd_i.W2_51[4]+=dW2_51.comps[4];
				bnd_i.W2_51[5]+=dW2_51.comps[5];
				bnd_i.W2_51[6]+=dW2_51.comps[6];
				W_rs_nu dW2_42=get_W2_rs_type_II_k(4,2,xk,xk_p1,nk,Lk,bnd_i.dt[j][k]);
				bnd_i.W2_42[0]+=dW2_42.comps[0];
				bnd_i.W2_42[1]+=dW2_42.comps[1];
				bnd_i.W2_42[2]+=dW2_42.comps[2];
				bnd_i.W2_42[3]+=dW2_42.comps[3];
				bnd_i.W2_42[4]+=dW2_42.comps[4];
				bnd_i.W2_42[5]+=dW2_42.comps[5];
				bnd_i.W2_42[6]+=dW2_42.comps[6];
				W_rs_nu dW2_33=get_W2_rs_type_II_k(3,3,xk,xk_p1,nk,Lk,bnd_i.dt[j][k]);
				bnd_i.W2_33[0]+=dW2_33.comps[0];
				bnd_i.W2_33[1]+=dW2_33.comps[1];
				bnd_i.W2_33[2]+=dW2_33.comps[2];
				bnd_i.W2_33[3]+=dW2_33.comps[3];
				bnd_i.W2_33[4]+=dW2_33.comps[4];
				bnd_i.W2_33[5]+=dW2_33.comps[5];
				bnd_i.W2_33[6]+=dW2_33.comps[6];
				W_rs_nu dW2_24=get_W2_rs_type_II_k(2,4,xk,xk_p1,nk,Lk,bnd_i.dt[j][k]);
				bnd_i.W2_24[0]+=dW2_24.comps[0];
				bnd_i.W2_24[1]+=dW2_24.comps[1];
				bnd_i.W2_24[2]+=dW2_24.comps[2];
				bnd_i.W2_24[3]+=dW2_24.comps[3];
				bnd_i.W2_24[4]+=dW2_24.comps[4];
				bnd_i.W2_24[5]+=dW2_24.comps[5];
				bnd_i.W2_24[6]+=dW2_24.comps[6];
				W_rs_nu dW2_15=get_W2_rs_type_II_k(1,5,xk,xk_p1,nk,Lk,bnd_i.dt[j][k]);
				bnd_i.W2_15[0]+=dW2_15.comps[0];
				bnd_i.W2_15[1]+=dW2_15.comps[1];
				bnd_i.W2_15[2]+=dW2_15.comps[2];
				bnd_i.W2_15[3]+=dW2_15.comps[3];
				bnd_i.W2_15[4]+=dW2_15.comps[4];
				bnd_i.W2_15[5]+=dW2_15.comps[5];
				bnd_i.W2_15[6]+=dW2_15.comps[6];
				W_rs_nu dW2_06=get_W2_rs_type_II_k(0,6,xk,xk_p1,nk,Lk,bnd_i.dt[j][k]);
				bnd_i.W2_06[0]+=dW2_06.comps[0];
				bnd_i.W2_06[1]+=dW2_06.comps[1];
				bnd_i.W2_06[2]+=dW2_06.comps[2];
				bnd_i.W2_06[3]+=dW2_06.comps[3];
				bnd_i.W2_06[4]+=dW2_06.comps[4];
				bnd_i.W2_06[5]+=dW2_06.comps[5];
				bnd_i.W2_06[6]+=dW2_06.comps[6];	
			}
		}
		//
		bnds.push_back(bnd_i);
	}
	return bnds;
}
// 
// write out appended bnd_solve coordinates
void write_bnd1(bnd bnd_i,string fname) {
	ofstream myfile;
	myfile.open(fname);
	for (size_t i=0;i<bnd_i.x.size();++i) {
		for (size_t j=0;j<bnd_i.x[i].size();++j) {
			myfile << bnd_i.x[i][j] << "," << bnd_i.y[i][j] << endl;
		}
	}
	myfile.close();
}
//
// write out non-appended bnd_solve coordinates
void write_bnd2(bnd bnd_i,string fname) {
	ofstream myfile;
	myfile.open(fname);
	for (size_t i=0;i<bnd_i.x.size();++i) {
		myfile << bnd_i.x[i].size();
		if (i<bnd_i.x.size()-1) {myfile << ",";}
	}
	myfile << endl;
	for (size_t i=0;i<bnd_i.x.size();++i) {
		for (size_t j=0;j<bnd_i.x[i].size();++j) {
			myfile << bnd_i.x[i][j] << "," << bnd_i.y[i][j] << endl;
		}
	}
	myfile.close();
}
//
// write out minkowski characterization results of type 0, i.e. w0_rs
void write_w0_rs(vector<bnd> bnds, string fname) {
	// 
	// exhaustive case of w0_rs up to order 6, i.e. r+s<=6
	// w0_00,
	// w0_10_x,w0_10_y,
	// w0_20_xx,w0_20_yy,w0_20_xy,
	// w0_30_xxx,w0_30_yyy,w0_30_xxy,w0_30_yyx,
	// w0_40_xxxx,w0_40_yyyy,w0_40_xxxy,w0_40_yyyx,w0_40_xxyy,
	// w0_50_xxxxx,w0_50_yyyyy,w0_50_xxxxy,w0_50_yyyyx,w0_50_xxxyy,w0_50_yyyxx,
	// w0_60_xxxxxx,w0_60_yyyyyy,w0_60_xxxxxy,w0_60_yyyyyx,w0_60_xxxxyy,w0_60_yyyyxx,w0_60_xxxyyy,
	//
	ofstream myfile;
	myfile.open(fname);
	myfile.precision(13);
	for (size_t i=0;i<bnds.size();++i) {
		myfile << bnds[i].W0 << ",";
		myfile << bnds[i].W0_10[0] << "," << bnds[i].W0_10[1] << ",";
		myfile << bnds[i].W0_20[0] << "," << bnds[i].W0_20[1] << "," << bnds[i].W0_20[2] << ",";
		myfile << bnds[i].W0_30[0] << "," << bnds[i].W0_30[1] << "," << bnds[i].W0_30[2] << "," << bnds[i].W0_30[3] << ",";
		myfile << bnds[i].W0_40[0] << "," << bnds[i].W0_40[1] << "," << bnds[i].W0_40[2] << "," << bnds[i].W0_40[3] << "," << bnds[i].W0_40[4] << ",";
		myfile << bnds[i].W0_50[0] << "," << bnds[i].W0_50[1] << "," << bnds[i].W0_50[2] << "," << bnds[i].W0_50[3] << "," << bnds[i].W0_50[4] << "," << bnds[i].W0_50[5] << ",";
		myfile << bnds[i].W0_60[0] << "," << bnds[i].W0_60[1] << "," << bnds[i].W0_60[2] << "," << bnds[i].W0_60[3] << "," << bnds[i].W0_60[4] << "," << bnds[i].W0_60[5] << "," << bnds[i].W0_60[6] << endl; 
	}
	myfile.close();
}
//
// write out minkowski characterization results of type 1, i.e. w1_rs
void write_w1_rs(vector<bnd> bnds, string fname) {
	// 
	// exhaustive case of w1_rs up to order 6, i.e. r+s<=6
	// w1_00,
	// w1_10_x,w1_10_y,
	// w1_20_xx,w1_20_yy,w1_20_xy,
	// w1_30_xxx,w1_30_yyy,w1_30_xxy,w1_30_yyx,
	// w1_40_xxxx,w1_40_yyyy,w1_40_xxxy,w1_40_yyyx,w1_40_xxyy,
	// w1_50_xxxxx,w1_50_yyyyy,w1_50_xxxxy,w1_50_yyyyx,w1_50_xxxyy,w1_50_yyyxx,
	// w1_60_xxxxxx,w1_60_yyyyyy,w1_60_xxxxxy,w1_60_yyyyyx,w1_60_xxxxyy,w1_60_yyyyxx,w1_60_xxxyyy,
	// w1_01_x,w1_01_y,
	// w1_02_xx,w1_02_yy,w1_02_xy,
	// w1_03_xxx,w1_03_yyy,w1_03_xxy,w1_03_yyx,
	// w1_04_xxxx,w1_04_yyyy,w1_04_xxxy,w1_04_yyyx,w1_04_xxyy,
	// w1_05_xxxxx,w1_05_yyyyy,w1_05_xxxxy,w1_05_yyyyx,w1_05_xxxyy,w1_05_yyyxx,
	// w1_06_xxxxxx,w1_06_yyyyyy,w1_06_xxxxxy,w1_06_yyyyyx,w1_06_xxxxyy,w1_06_yyyyxx,w1_06_xxxyyy,
	// w1_11_xx,w1_11_yy,w1_11_xy,
	// w1_21_xxx,w1_21_yyy,w1_21_xxy,w1_21_yyx,
	// w1_12_xxx,w1_12_yyy,w1_12_xxy,w1_12_yyx,
	// w1_31_xxxx,w1_31_yyyy,w1_31_xxxy,w1_31_yyyx,w1_31_xxyy,
	// w1_22_xxxx,w1_22_yyyy,w1_22_xxxy,w1_22_yyyx,w1_22_xxyy,
	// w1_13_xxxx,w1_13_yyyy,w1_13_xxxy,w1_13_yyyx,w1_13_xxyy.
	// w1_41_xxxxx,w1_41_yyyyy,w1_41_xxxxy,w1_41_yyyyx,w1_41_xxxyy,w1_41_yyyxx,
	// w1_32_xxxxx,w1_32_yyyyy,w1_32_xxxxy,w1_32_yyyyx,w1_32_xxxyy,w1_32_yyyxx,
	// w1_23_xxxxx,w1_23_yyyyy,w1_23_xxxxy,w1_23_yyyyx,w1_23_xxxyy,w1_23_yyyxx,
	// w1_14_xxxxx,w1_14_yyyyy,w1_14_xxxxy,w1_14_yyyyx,w1_14_xxxyy,w1_14_yyyxx,
	// w1_51_xxxxxx,w1_51_yyyyyy,w1_51_xxxxxy,w1_51_yyyyyx,w1_51_xxxxyy,w1_51_yyyyxx,w1_51_xxxyyy,
	// w1_42_xxxxxx,w1_42_yyyyyy,w1_42_xxxxxy,w1_42_yyyyyx,w1_42_xxxxyy,w1_42_yyyyxx,w1_42_xxxyyy,
	// w1_33_xxxxxx,w1_33_yyyyyy,w1_33_xxxxxy,w1_33_yyyyyx,w1_33_xxxxyy,w1_33_yyyyxx,w1_33_xxxyyy,
	// w1_24_xxxxxx,w1_24_yyyyyy,w1_24_xxxxxy,w1_24_yyyyyx,w1_24_xxxxyy,w1_24_yyyyxx,w1_24_xxxyyy,
	// w1_15_xxxxxx,w1_15_yyyyyy,w1_15_xxxxxy,w1_15_yyyyyx,w1_15_xxxxyy,w1_15_yyyyxx,w1_15_xxxyyy.
	//
	ofstream myfile;
	myfile.open(fname);
	myfile.precision(13);
	for (size_t i=0;i<bnds.size();++i) {
		myfile << bnds[i].W1 << ",";
		myfile << bnds[i].W1_10[0] << "," << bnds[i].W1_10[1] << ",";
		myfile << bnds[i].W1_20[0] << "," << bnds[i].W1_20[1] << "," << bnds[i].W1_20[2] << ",";
		myfile << bnds[i].W1_30[0] << "," << bnds[i].W1_30[1] << "," << bnds[i].W1_30[2] << "," << bnds[i].W1_30[3] << ",";
		myfile << bnds[i].W1_40[0] << "," << bnds[i].W1_40[1] << "," << bnds[i].W1_40[2] << "," << bnds[i].W1_40[3] << "," << bnds[i].W1_40[4] << ",";
		myfile << bnds[i].W1_50[0] << "," << bnds[i].W1_50[1] << "," << bnds[i].W1_50[2] << "," << bnds[i].W1_50[3] << "," << bnds[i].W1_50[4] << "," << bnds[i].W1_50[5] << ",";
		myfile << bnds[i].W1_60[0] << "," << bnds[i].W1_60[1] << "," << bnds[i].W1_60[2] << "," << bnds[i].W1_60[3] << "," << bnds[i].W1_60[4] << "," << bnds[i].W1_60[5] << "," << bnds[i].W1_60[6] << ",";
		myfile << bnds[i].W1_01[0] << "," << bnds[i].W1_01[1] << ",";
		myfile << bnds[i].W1_02[0] << "," << bnds[i].W1_02[1] << "," << bnds[i].W1_02[2] << ",";
		myfile << bnds[i].W1_03[0] << "," << bnds[i].W1_03[1] << "," << bnds[i].W1_03[2] << "," << bnds[i].W1_03[3] << ",";
		myfile << bnds[i].W1_04[0] << "," << bnds[i].W1_04[1] << "," << bnds[i].W1_04[2] << "," << bnds[i].W1_04[3] << "," << bnds[i].W1_04[4] << ",";
		myfile << bnds[i].W1_05[0] << "," << bnds[i].W1_05[1] << "," << bnds[i].W1_05[2] << "," << bnds[i].W1_05[3] << "," << bnds[i].W1_05[4] << "," << bnds[i].W1_05[5] << ",";
		myfile << bnds[i].W1_06[0] << "," << bnds[i].W1_06[1] << "," << bnds[i].W1_06[2] << "," << bnds[i].W1_06[3] << "," << bnds[i].W1_06[4] << "," << bnds[i].W1_06[5] << "," << bnds[i].W1_06[6] << ",";
		myfile << bnds[i].W1_11[0] << "," << bnds[i].W1_11[1] << "," << bnds[i].W1_11[2] << ",";
		myfile << bnds[i].W1_21[0] << "," << bnds[i].W1_21[1] << "," << bnds[i].W1_21[2] << "," << bnds[i].W1_21[3] << ",";
		myfile << bnds[i].W1_12[0] << "," << bnds[i].W1_12[1] << "," << bnds[i].W1_12[2] << "," << bnds[i].W1_12[3] << ",";
		myfile << bnds[i].W1_31[0] << "," << bnds[i].W1_31[1] << "," << bnds[i].W1_31[2] << "," << bnds[i].W1_31[3] << "," << bnds[i].W1_31[4] << ",";
		myfile << bnds[i].W1_22[0] << "," << bnds[i].W1_22[1] << "," << bnds[i].W1_22[2] << "," << bnds[i].W1_22[3] << "," << bnds[i].W1_22[4] << ",";
		myfile << bnds[i].W1_13[0] << "," << bnds[i].W1_13[1] << "," << bnds[i].W1_13[2] << "," << bnds[i].W1_13[3] << "," << bnds[i].W1_13[4] << ",";
		myfile << bnds[i].W1_41[0] << "," << bnds[i].W1_41[1] << "," << bnds[i].W1_41[2] << "," << bnds[i].W1_41[3] << "," << bnds[i].W1_41[4] << "," << bnds[i].W1_41[5] << ",";
		myfile << bnds[i].W1_32[0] << "," << bnds[i].W1_32[1] << "," << bnds[i].W1_32[2] << "," << bnds[i].W1_32[3] << "," << bnds[i].W1_32[4] << "," << bnds[i].W1_32[5] << ",";
		myfile << bnds[i].W1_23[0] << "," << bnds[i].W1_23[1] << "," << bnds[i].W1_23[2] << "," << bnds[i].W1_23[3] << "," << bnds[i].W1_23[4] << "," << bnds[i].W1_23[5] << ",";
		myfile << bnds[i].W1_14[0] << "," << bnds[i].W1_14[1] << "," << bnds[i].W1_14[2] << "," << bnds[i].W1_14[3] << "," << bnds[i].W1_14[4] << "," << bnds[i].W1_14[5] << ",";
		myfile << bnds[i].W1_51[0] << "," << bnds[i].W1_51[1] << "," << bnds[i].W1_51[2] << "," << bnds[i].W1_51[3] << "," << bnds[i].W1_51[4] << "," << bnds[i].W1_51[5] << "," << bnds[i].W1_51[6] << ",";
		myfile << bnds[i].W1_42[0] << "," << bnds[i].W1_42[1] << "," << bnds[i].W1_42[2] << "," << bnds[i].W1_42[3] << "," << bnds[i].W1_42[4] << "," << bnds[i].W1_42[5] << "," << bnds[i].W1_42[6] << ",";
		myfile << bnds[i].W1_33[0] << "," << bnds[i].W1_33[1] << "," << bnds[i].W1_33[2] << "," << bnds[i].W1_33[3] << "," << bnds[i].W1_33[4] << "," << bnds[i].W1_33[5] << "," << bnds[i].W1_33[6] << ",";
		myfile << bnds[i].W1_24[0] << "," << bnds[i].W1_24[1] << "," << bnds[i].W1_24[2] << "," << bnds[i].W1_24[3] << "," << bnds[i].W1_24[4] << "," << bnds[i].W1_24[5] << "," << bnds[i].W1_24[6] << ",";
		myfile << bnds[i].W1_15[0] << "," << bnds[i].W1_15[1] << "," << bnds[i].W1_15[2] << "," << bnds[i].W1_15[3] << "," << bnds[i].W1_15[4] << "," << bnds[i].W1_15[5] << "," << bnds[i].W1_15[6] << endl; 
	}
	myfile.close();
}
//
// write out minkowski characterization results of type 1, i.e. w2_rs
void write_w2_rs(vector<bnd> bnds, string fname) {
	// 
	// exhaustive case of W2_rs up to order 6, i.e. r+s<=6
	// W2_00,
	// W2_10_x,W2_10_y,
	// W2_20_xx,W2_20_yy,W2_20_xy,
	// W2_30_xxx,W2_30_yyy,W2_30_xxy,W2_30_yyx,
	// W2_40_xxxx,W2_40_yyyy,W2_40_xxxy,W2_40_yyyx,W2_40_xxyy,
	// W2_50_xxxxx,W2_50_yyyyy,W2_50_xxxxy,W2_50_yyyyx,W2_50_xxxyy,W2_50_yyyxx,
	// W2_60_xxxxxx,W2_60_yyyyyy,W2_60_xxxxxy,W2_60_yyyyyx,W2_60_xxxxyy,W2_60_yyyyxx,W2_60_xxxyyy,
	// W2_01_x,W2_01_y,
	// W2_02_xx,W2_02_yy,W2_02_xy,
	// W2_03_xxx,W2_03_yyy,W2_03_xxy,W2_03_yyx,
	// W2_04_xxxx,W2_04_yyyy,W2_04_xxxy,W2_04_yyyx,W2_04_xxyy,
	// W2_05_xxxxx,W2_05_yyyyy,W2_05_xxxxy,W2_05_yyyyx,W2_05_xxxyy,W2_05_yyyxx,
	// W2_06_xxxxxx,W2_06_yyyyyy,W2_06_xxxxxy,W2_06_yyyyyx,W2_06_xxxxyy,W2_06_yyyyxx,W2_06_xxxyyy,
	// W2_11_xx,W2_11_yy,W2_11_xy,
	// W2_21_xxx,W2_21_yyy,W2_21_xxy,W2_21_yyx,
	// W2_12_xxx,W2_12_yyy,W2_12_xxy,W2_12_yyx,
	// W2_31_xxxx,W2_31_yyyy,W2_31_xxxy,W2_31_yyyx,W2_31_xxyy,
	// W2_22_xxxx,W2_22_yyyy,W2_22_xxxy,W2_22_yyyx,W2_22_xxyy,
	// W2_13_xxxx,W2_13_yyyy,W2_13_xxxy,W2_13_yyyx,W2_13_xxyy.
	// W2_41_xxxxx,W2_41_yyyyy,W2_41_xxxxy,W2_41_yyyyx,W2_41_xxxyy,W2_41_yyyxx,
	// W2_32_xxxxx,W2_32_yyyyy,W2_32_xxxxy,W2_32_yyyyx,W2_32_xxxyy,W2_32_yyyxx,
	// W2_23_xxxxx,W2_23_yyyyy,W2_23_xxxxy,W2_23_yyyyx,W2_23_xxxyy,W2_23_yyyxx,
	// W2_14_xxxxx,W2_14_yyyyy,W2_14_xxxxy,W2_14_yyyyx,W2_14_xxxyy,W2_14_yyyxx,
	// W2_51_xxxxxx,W2_51_yyyyyy,W2_51_xxxxxy,W2_51_yyyyyx,W2_51_xxxxyy,W2_51_yyyyxx,W2_51_xxxyyy,
	// W2_42_xxxxxx,W2_42_yyyyyy,W2_42_xxxxxy,W2_42_yyyyyx,W2_42_xxxxyy,W2_42_yyyyxx,W2_42_xxxyyy,
	// W2_33_xxxxxx,W2_33_yyyyyy,W2_33_xxxxxy,W2_33_yyyyyx,W2_33_xxxxyy,W2_33_yyyyxx,W2_33_xxxyyy,
	// W2_24_xxxxxx,W2_24_yyyyyy,W2_24_xxxxxy,W2_24_yyyyyx,W2_24_xxxxyy,W2_24_yyyyxx,W2_24_xxxyyy,
	// W2_15_xxxxxx,W2_15_yyyyyy,W2_15_xxxxxy,W2_15_yyyyyx,W2_15_xxxxyy,W2_15_yyyyxx,W2_15_xxxyyy.
	//
	ofstream myfile;
	myfile.open(fname);
	myfile.precision(13);
	for (size_t i=0;i<bnds.size();++i) {
		myfile << bnds[i].W2 << ",";
		myfile << bnds[i].W2_10[0] << "," << bnds[i].W2_10[1] << ",";
		myfile << bnds[i].W2_20[0] << "," << bnds[i].W2_20[1] << "," << bnds[i].W2_20[2] << ",";
		myfile << bnds[i].W2_30[0] << "," << bnds[i].W2_30[1] << "," << bnds[i].W2_30[2] << "," << bnds[i].W2_30[3] << ",";
		myfile << bnds[i].W2_40[0] << "," << bnds[i].W2_40[1] << "," << bnds[i].W2_40[2] << "," << bnds[i].W2_40[3] << "," << bnds[i].W2_40[4] << ",";
		myfile << bnds[i].W2_50[0] << "," << bnds[i].W2_50[1] << "," << bnds[i].W2_50[2] << "," << bnds[i].W2_50[3] << "," << bnds[i].W2_50[4] << "," << bnds[i].W2_50[5] << ",";
		myfile << bnds[i].W2_60[0] << "," << bnds[i].W2_60[1] << "," << bnds[i].W2_60[2] << "," << bnds[i].W2_60[3] << "," << bnds[i].W2_60[4] << "," << bnds[i].W2_60[5] << "," << bnds[i].W2_60[6] << ",";
		myfile << bnds[i].W2_01[0] << "," << bnds[i].W2_01[1] << ",";
		myfile << bnds[i].W2_02[0] << "," << bnds[i].W2_02[1] << "," << bnds[i].W2_02[2] << ",";
		myfile << bnds[i].W2_03[0] << "," << bnds[i].W2_03[1] << "," << bnds[i].W2_03[2] << "," << bnds[i].W2_03[3] << ",";
		myfile << bnds[i].W2_04[0] << "," << bnds[i].W2_04[1] << "," << bnds[i].W2_04[2] << "," << bnds[i].W2_04[3] << "," << bnds[i].W2_04[4] << ",";
		myfile << bnds[i].W2_05[0] << "," << bnds[i].W2_05[1] << "," << bnds[i].W2_05[2] << "," << bnds[i].W2_05[3] << "," << bnds[i].W2_05[4] << "," << bnds[i].W2_05[5] << ",";
		myfile << bnds[i].W2_06[0] << "," << bnds[i].W2_06[1] << "," << bnds[i].W2_06[2] << "," << bnds[i].W2_06[3] << "," << bnds[i].W2_06[4] << "," << bnds[i].W2_06[5] << "," << bnds[i].W2_06[6] << ",";
		myfile << bnds[i].W2_11[0] << "," << bnds[i].W2_11[1] << "," << bnds[i].W2_11[2] << ",";
		myfile << bnds[i].W2_21[0] << "," << bnds[i].W2_21[1] << "," << bnds[i].W2_21[2] << "," << bnds[i].W2_21[3] << ",";
		myfile << bnds[i].W2_12[0] << "," << bnds[i].W2_12[1] << "," << bnds[i].W2_12[2] << "," << bnds[i].W2_12[3] << ",";
		myfile << bnds[i].W2_31[0] << "," << bnds[i].W2_31[1] << "," << bnds[i].W2_31[2] << "," << bnds[i].W2_31[3] << "," << bnds[i].W2_31[4] << ",";
		myfile << bnds[i].W2_22[0] << "," << bnds[i].W2_22[1] << "," << bnds[i].W2_22[2] << "," << bnds[i].W2_22[3] << "," << bnds[i].W2_22[4] << ",";
		myfile << bnds[i].W2_13[0] << "," << bnds[i].W2_13[1] << "," << bnds[i].W2_13[2] << "," << bnds[i].W2_13[3] << "," << bnds[i].W2_13[4] << ",";
		myfile << bnds[i].W2_41[0] << "," << bnds[i].W2_41[1] << "," << bnds[i].W2_41[2] << "," << bnds[i].W2_41[3] << "," << bnds[i].W2_41[4] << "," << bnds[i].W2_41[5] << ",";
		myfile << bnds[i].W2_32[0] << "," << bnds[i].W2_32[1] << "," << bnds[i].W2_32[2] << "," << bnds[i].W2_32[3] << "," << bnds[i].W2_32[4] << "," << bnds[i].W2_32[5] << ",";
		myfile << bnds[i].W2_23[0] << "," << bnds[i].W2_23[1] << "," << bnds[i].W2_23[2] << "," << bnds[i].W2_23[3] << "," << bnds[i].W2_23[4] << "," << bnds[i].W2_23[5] << ",";
		myfile << bnds[i].W2_14[0] << "," << bnds[i].W2_14[1] << "," << bnds[i].W2_14[2] << "," << bnds[i].W2_14[3] << "," << bnds[i].W2_14[4] << "," << bnds[i].W2_14[5] << ",";
		myfile << bnds[i].W2_51[0] << "," << bnds[i].W2_51[1] << "," << bnds[i].W2_51[2] << "," << bnds[i].W2_51[3] << "," << bnds[i].W2_51[4] << "," << bnds[i].W2_51[5] << "," << bnds[i].W2_51[6] << ",";
		myfile << bnds[i].W2_42[0] << "," << bnds[i].W2_42[1] << "," << bnds[i].W2_42[2] << "," << bnds[i].W2_42[3] << "," << bnds[i].W2_42[4] << "," << bnds[i].W2_42[5] << "," << bnds[i].W2_42[6] << ",";
		myfile << bnds[i].W2_33[0] << "," << bnds[i].W2_33[1] << "," << bnds[i].W2_33[2] << "," << bnds[i].W2_33[3] << "," << bnds[i].W2_33[4] << "," << bnds[i].W2_33[5] << "," << bnds[i].W2_33[6] << ",";
		myfile << bnds[i].W2_24[0] << "," << bnds[i].W2_24[1] << "," << bnds[i].W2_24[2] << "," << bnds[i].W2_24[3] << "," << bnds[i].W2_24[4] << "," << bnds[i].W2_24[5] << "," << bnds[i].W2_24[6] << ",";
		myfile << bnds[i].W2_15[0] << "," << bnds[i].W2_15[1] << "," << bnds[i].W2_15[2] << "," << bnds[i].W2_15[3] << "," << bnds[i].W2_15[4] << "," << bnds[i].W2_15[5] << "," << bnds[i].W2_15[6] << endl; 
	}
	myfile.close();
}
//
// testing for inf and nan values
void testing(vector<bnd> bnds) {
	for (size_t i=0;i<bnds.size();++i) {
		for (size_t j=0;j<bnds[i].x.size();++j) {
			for (size_t k=0;k<bnds[i].x[j].size();++k) {
				if (isinf(bnds[i].x[j][k])||isnan(bnds[i].x[j][k])) {
					cout << "i= " << i << ", j= " << j << endl;
				}
				else if (isinf(bnds[i].y[j][k])||isnan(bnds[i].y[j][k])) {
					cout << "i= " << i << ", j= " << j << endl;
				}
			}
		}
	}	
}

////////////////////////////////  main
// 
int main(int argc, char* argv[] ) {	
	string fname_in1=argv[1]; // file name for bnd_solve_data
	string fname_in2=argv[2]; // file name for mpp
	string fname_in3=argv[3]; // file name for neighbs
	int npts=atoi(argv[4]); // number of points per segment of common curve
	//
	vector<int> badgrain_ids={35,48,50,62,88,94,99,59,72,78}; // igrain_ids of cliques not solvable by bnd_solve
	//
	bnd_solve_data data=import_bnd_solve_data(fname_in1); // read data from bnd_solve
	marked_point_pattern egs_mpp=read_mpp(fname_in2); // read data from marked point pattern
	vector<vector<int>> neighbs=read_neighbs(fname_in3); // read data from list of cliques from gtess
	vector<clique> cliques=create_cliques(data,egs_mpp,neighbs); // create list of all cliques solvable by bnd_solve
	//
	vector<bnd> bnds=resolve_bnds(cliques,npts); // resolve for bnds of each solvable cliques
	testing(bnds); // testing for inf and nan values
	//
	int bnd_out_flag=0;
	// write bnd_resolve out
	if (bnd_out_flag==1) {
		//cout.precision(10);
		for (size_t i=0;i<bnds.size();++i) {
			ostringstream num;
			num << i;
			write_bnd1(bnds[i],"bnd_cpp_"+num.str()+".out"); // append common curves
		}
	}	
	//
	ostringstream num;
	num << npts;	
	write_w1_rs(bnds,"w1_rs_"+num.str()+".out");
	write_w2_rs(bnds,"w2_rs_"+num.str()+".out");
	write_w0_rs(bnds,"w0_rs_"+num.str()+".out");
	//
	return 0;
}

// Typical call:
// ./bnd_resolve bnd_solve_data.out mppsim_test01.out gtess_egt_test01.out 6000 






