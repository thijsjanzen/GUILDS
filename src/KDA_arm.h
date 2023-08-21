#ifndef KDA_ARM_H
#define KDA_ARM_H

#include <cmath>
#include <vector>
#include <algorithm>

class log_val {
public:
  log_val() {
    is_zero_ = true;
    log_val_ = -1e10;
  }

  log_val(double v) : log_val_(v) {
    is_zero_ = false;
  }

  void reset() {
    is_zero_ = true;
    log_val_ = -1e10;
  }

  void set_with_norm(double v) {
    log_val_ = log(v);
    is_zero_ = false;
  }

  void set_with_log(double l_v) {
    log_val_ = l_v;
    is_zero_ = false;
  }

  bool is_zero() const {
    return is_zero_;
  }

  double get_log_val() const {
    return log_val_;
  }

  log_val& operator+=(const log_val& other) {
    if (is_zero_) {
      if (other.is_zero()) {
        return *this;
      }
      log_val_ = other.get_log_val();
      is_zero_ = false;
      return *this;
    }


    auto l_a = this->get_log_val();
    auto l_b = other.get_log_val();
    log_val_ = l_a + log(1.0 + exp(l_b - l_a));
    return *this;
  }

private:
  double log_val_;
  bool is_zero_;
};

log_val operator+(const log_val& a, const log_val& b) {
  if (a.is_zero()) {
    return b;
  }
  if (b.is_zero()) {
    return a;
  }
  auto l_a = a.get_log_val();
  auto l_b = b.get_log_val();
  auto res = l_a + log(1.0 + exp(l_b - l_a));
  return log_val(res);
}

double log(const log_val& a) {
  if (a.is_zero()) {
    return log(0.0);
  }
  return a.get_log_val();
}

log_val operator*(const log_val& a, const log_val& b) {
  double l_a = a.get_log_val();
  double l_b = b.get_log_val();
  double res = l_a + l_b;
  return log_val(res);
}

std::vector<double> calcLogKDA_arm(
    size_t numspecies,
    std::vector<size_t> Abund)  {


  if (Abund.size() < 1) return std::vector<double>();

  std::sort(Abund.begin(), Abund.end());

  double J = 0;
  double SPP = numspecies;
  for(size_t s = 0; s < SPP; s++) {
    J += Abund[s];
  }

  double MaxA = 0;
  MaxA = Abund.back(); //MaxA = Abund[(int)SPP-1];

  // abundance distribution
  std::vector< size_t > Phi(MaxA + 1, 0); //int *Phi = new int[MaxA+1]; for(s=0;s<=MaxA;s++) Phi[s]=0;

  for (size_t s = 0; s < SPP; s++) Phi[ Abund[s] ]++;

  // Number of distinct abundances
  size_t NDA = 0;
  for (size_t s = 0; s <= MaxA; s++) if(Phi[s] > 0) {NDA++;}


  //cerr << "Start computing Stirling numbers ...\n";
  // FIRST STAGE: compute the Stirling numbers
  // I use the relationship S(n,m)=S(n-1,m-1)-(n-1)S(n-1,m)
  // where n >= m >= 1
  // The relation of interest is sum_m S(n,m)/S(n,1)*S(m,1) * x^m
  // defined to be equal to sum_m T(n,m) * x^m
  // The recurrence relation on T(n,m) is
  // T(n,m)= T(n-1,m) + T(n-1,m-1)*(m-1)/(n-1)

  std::vector<size_t> f(NDA, 0); //int *f = new int[NDA];
  std::vector<size_t> g(NDA, 0); //int *g = new int[NDA];
  size_t i = 0;
  //for(s=0;s<NDA;s++) {f[s]=0;g[s]=0;}
  for(int n=0; n <= MaxA; n++) {
    if(Phi[n] > 0) {
      f[i] = Phi[n];
      g[i] = n;
      i++;
    }
  }

  std::vector< double > fill;
  std::vector< std::vector< double > > T(NDA,fill);

  std::vector< double > T0(g[0]+1);
  T[0]    = T0;
  T[0][0] = 0;
  T[0][1] = 1;

  if(g[0]!=1) {
    std::vector<double> lS2(g[0]+1); //        double *lS2 = new double[g[0]+1];
    lS2[0]=0;lS2[1]=1;
    for (size_t n=2;n<=g[0];n++) {
      std::vector<double> lS1(n+1); //double *lS1 = new double[n+1];
      for(size_t im=0;im<=n-1;im++) {
        lS1[im] = lS2[im];
      }
      lS1[n]=0;
      for(size_t im=2;im<=n;im++) {
        lS2[im] = lS1[im]+lS1[im-1]*(im-1)/(n-1);
      }
    }
    for(size_t im=2;im<=g[0];im++) {
      T[0][im]=lS2[im];
    }
  }


  for (size_t in = 1; in < i; in++) {
    std::vector<double> Tin(g[in]+1); //T[in]= new double[g[in]+1];
    T[in] = Tin;
    T[in][0]=0;
    T[in][1]=1;

    std::vector<double> lS2(g[in]+1); //double *lS2 = new double[g[in]+1];
    for(size_t im=0;im<=g[in-1];im++) {
      lS2[im] = T[in-1][im];
    }
    for (size_t n=g[in-1]+1;n<=g[in];n++) {
      std::vector<double> lS1(n+1); //double *lS1 = new double[n+1];
      for(size_t im=0;im<=n-1;im++) {
        lS1[im] = lS2[im];
      }
      lS1[n]=0;
      for(size_t im=2;im<=n;im++) {
        lS2[im] = lS1[im]+lS1[im-1]*(im-1)/(n-1);
      }
    }
    for(size_t im=2;im<=g[in];im++) {
      T[in][im]=lS2[im];
    }
  }
  // After this stage we have stored in T[i][m] T(g[i],m)
  // with T(n,m) = S(n,m)*S(m,1)/S(n,1) for i>0

  // precompute log values
  for (auto& i : T) {
    for (auto& j : i) {
      j = log(j);
    }
  }

  // cerr << "Start computing ln(K(D,A)) ...\n";
  // SECOND STAGE: compute the K(D,A)
  // I follow Etienne's route. Compute the product of polynomials
  // of length J

  std::vector< log_val > K_prime(J+2);
  std::vector< log_val > log_poly2(J+1);

  K_prime[0].set_with_norm(1);

  size_t degree = 0;
  size_t spe = 0;
  double minus = 4500.0 * log(10) / SPP;

  for (size_t i=0;i<NDA;i++) { // loop over number of distinct abundances
    for (size_t j=0;j<f[i];j++) { // loop over abundances per class

      for (size_t nn=0; nn<= degree; nn++) {
        for (size_t mm = 1; mm <= g[i]; mm++){
          if (!K_prime[nn].is_zero()){
            log_poly2[nn + mm] += T[i][mm] + K_prime[nn].get_log_val();
          }
        }
      }

      degree += g[i];
      for(size_t nn=0;nn<=degree;nn++){
        if (!log_poly2[nn].is_zero() ) {
          double add_new = log_poly2[nn].get_log_val() - minus;
          K_prime[nn].set_with_log(add_new);
        } else {
          K_prime[nn].reset();
        }

        log_poly2[nn].reset();
      }
      spe++;
    }
  }


  std::vector<double> K(J+2, 0.0);

  size_t cnt = 0;
  for (const auto& i : K_prime) {
    K[cnt] = static_cast<double>(i.get_log_val());
    if (K[cnt] == -1e10) K[cnt] = 0.0;
    cnt++;
  }

  double borneinf = (SPP - 1);
  double bornesup = (J + 1 );
  double maxlog = 11333.2;
  int infinity = 0;
  for(int i= SPP; i <= J; i++){
    if (i > (int)K.size()) {
      return K; //throw std::out_of_range("i > K.size()");
    }

    if ((K[i]>maxlog) || (K[i] < -maxlog)) {
      infinity = 1;
      break;
    }
    borneinf++;
  }  //after that, borneinf=indice next to infinity but before

  for(int i = 0; i <= J - SPP; i++){
    if (i > (int)K.size()) {
      //throw std::out_of_range("i > K.size()");
      return K;
    }

    if ((K[J-i]>maxlog)||(K[(int)J-i]<-maxlog)) {
      infinity=1;
      break;
    }
    bornesup--;
  }    //after that, bornesup=indice next to infinity but after
  if ( infinity==1 ) {
    //cerr << "WARNING : the sample is too large to compute an exact likelihood, the program is thus doing approximations. The following results are to be taken with caution"<<endl;
    //cerr << "Value of A above which K(D,A) is computed approximately ="<<borneinf<<endl;
    //cerr << "Value of A below which K(D,A) is computed approximately ="<<bornesup<<endl;

    //fitting of the infinite values of K[A] by a polynom of degree 3
    //computing of the derivatives at the critic points

    if((int)borneinf >= (int)K.size()) {
      return K; //throw std::out_of_range("borneinf > K.size()");
    }
    if(bornesup >= (int)(K.size())) {
      return K; //throw std::out_of_range("bornesup > K.size()");
    }

    int index1 = (int)borneinf;
    int index2 = (int)borneinf - 1;

    if (index1 >= (int)K.size() || index1 < 0) {
      return K; //throw std::out_of_range("borneinf");
    }
    if (index2 >= (int)K.size() || index2 < 0) {
      return K; //throw std::out_of_range("borneinf - 1");
    }

    double Kprimeinf = K[index1] - K[index2];


    double Kprimesup1 = -50.0;  //out of range check

    if(((int)bornesup+1) < (int)K.size() &&
       bornesup+1 >= 0) Kprimesup1 = K[bornesup+1];

    double Kprimesup2 = 0.0;   //out of range check

    if((int)bornesup >= 0.0) Kprimesup2 = K[bornesup];

    double Kprimesup = Kprimesup1 - Kprimesup2; //K[bornesup+1]-K[bornesup];
    // definition of the parameters of the fitted polynom aX^3+bX^2+cX+d
    double a,b,c,d;
    //inversion of the linear system of equations (with the Gauss method)
    double borneinf2=(double)borneinf*(double)borneinf;
    double borneinf3=(double)borneinf2*(double)borneinf;
    double bornesup2=(double)bornesup*(double)bornesup;
    double bornesup3=(double)bornesup2*(double)bornesup;
    d=(Kprimesup - 3*bornesup2 * K[borneinf]/borneinf3+(2*bornesup/(double)borneinf-3*bornesup2/borneinf2)*(Kprimeinf-3*K[borneinf]/(double)borneinf)-((1+3*bornesup2/borneinf2-4*bornesup/(double)borneinf)/(bornesup-2*bornesup2/(double)borneinf+bornesup3/borneinf2))*(K[bornesup]-bornesup3*K[borneinf]/borneinf3+(bornesup2/(double)borneinf-bornesup3/borneinf2)*(Kprimeinf-3*K[borneinf]/(double)borneinf)))/((6*bornesup2/borneinf3)-(6*bornesup/borneinf2)-((1+3*bornesup2/borneinf2-4*bornesup/(double)borneinf)/(bornesup-2*bornesup2/(double)borneinf+bornesup3/borneinf2))*(1-3*bornesup2/borneinf2+2*bornesup3/borneinf3));
    c=((K[bornesup]-bornesup3 * K[borneinf]/borneinf3+(bornesup2/(double)borneinf-bornesup3/borneinf2)*(Kprimeinf-3*K[borneinf]/(double)borneinf))-d*(1-3*bornesup2/borneinf2+2*bornesup3/borneinf3))/(bornesup-2*bornesup2/(double)borneinf+bornesup3/borneinf2);
    b =
      (Kprimeinf -
      3*K[borneinf] /
        (double)borneinf +
          2*c +
          3*d / (double)borneinf) / (0.0 - (double)borneinf);
    a= (K[borneinf] -
      b*borneinf2 -
      c*(double)borneinf - d) /
        borneinf3;

    //reconstruction of K[A] with the fitted polynom
    for (int i=borneinf+1; i<bornesup; i++) {
      K[i]=(a*i*i*i+b*i*i+c*i+d);
    }
  }

  return K;
}



#endif
