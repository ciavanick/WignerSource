class wignerfunction
{
private:
    float m_hcut = 197.3; //MeV * fem
    int m_ndim; //number of dimensions
    float m_minX, m_maxX, m_minP, m_maxP; //lower and upper range limits
    float m_dx=0.1,m_dp=1; 
    float m_rmin, m_rmax, m_pmin, m_pmax; //range of action
    float m_r0; 
    float m_k;


    double m_const; //normalization costant

    TString m_name;

    TF2 *m_w = nullptr; //wigner function

    // wigner function definition
    static double wignerSource(double *x, double *pm)
    {
        double r = x[0];
        double p = x[1];

        return pm[0] * TMath::Exp(-r * r * 0.5 / (pm[1] * pm[1]) - ((p - pm[2]) * (p - pm[2])) * 0.5 * (pm[1] * pm[1]) / (pm[3] * pm[3])) / TMath::Power(2 * TMath::Pi() * pm[3], pm[4]);
    }

    // wigner function*jacobian
    static double JacobianFun(double *x,double *pm){
        double r = x[0];
        double p = x[1];
    
        float jacobian=TMath::Power(r*p,pm[4]-1);
        if(pm[4]==2) jacobian *= 2*TMath::Pi()*2*TMath::Pi();
        if(pm[4]==3) jacobian *= 4*TMath::Pi()*2*TMath::Pi() * (1-TMath::Exp(-2*pm[1]*pm[1]/(pm[3]*pm[3])*p*pm[2]))/p/pm[2]/(pm[1]*pm[1])*pm[3]*pm[3];
    
        return wignerSource(x,pm) * jacobian;
    }

    static double Jacobianw2(double *x,double *pm){
        double r = x[0];
        double p = x[1];
    
        float jacobian=TMath::Power(r*p,pm[4]-1);
        if(pm[4]==2) jacobian *= 2*TMath::Pi()*2*TMath::Pi();
        if(pm[4]==3) jacobian *= 4*TMath::Pi()*2*TMath::Pi() * (1-TMath::Exp(-4*pm[1]*pm[1]/(pm[3]*pm[3])*p*pm[2]))/p/pm[2]/(pm[1]*pm[1])*pm[3]*pm[3]*0.5;
    
        return wignerSource(x,pm) * jacobian;
    }

    //w * w with the jacobian included
    static double wignerSource2(double *x,double *pm){
        double r = x[0];
        double p = x[1];
    
        return wignerSource(x,pm) * Jacobianw2(x,pm);
    }

     //integration method
     double myintegral(TF2 *f){
        double res = 0;
        for(float x=m_dx/2+m_minX; x < m_maxX;x+=m_dx){
          for(float p=m_dp/2+m_minP; p < m_maxP;p+=m_dp){
            res += f->Eval(x,p);
          }
        }
        res *= m_dx*m_dp;
        return res;
      }

      double normalization(){
        TString jacobname = m_name + "jacob";

        //to compute the normalization a new TF2 is created and is just W x jacobian, then is integrated with myintegral
        TF2* jacobw1 = new TF2(jacobname, JacobianFun,  m_minX, m_maxX, m_minP, m_maxP, 5);
        jacobw1->SetParameter(0, 1.);
        jacobw1->SetParameter(1, m_r0);
        jacobw1->SetParameter(2, m_k);
        jacobw1->SetParameter(3, m_hcut);
        jacobw1->SetParameter(4, m_ndim);

        return 1./myintegral(jacobw1);
        
      }

public:
    //constructor definition
    wignerfunction(float r0, float k, int dim, float rmin, float rmax, float pmin, float pmax, TString name){
        m_ndim = dim;
        m_minX = -100*(m_ndim==1);
        m_maxX = 100;
        m_minP=-1000*(m_ndim==1);
        m_maxP=1000;
        m_name = name;

        m_r0 = r0;
        m_k = k;

        //defining the wigner function
        m_w = new TF2(m_name, wignerSource, rmin, rmax, pmin, pmax, 5);
        m_w->SetParameter(0, 1.);
        m_w->SetParameter(1, r0);
        m_w->SetParameter(2, k);
        m_w->SetParameter(3, m_hcut);
        m_w->SetParameter(4, m_ndim);

        m_const = normalization();
        m_w->SetParameter(0, m_const);

    }
    //method to draw W
    void const DrawWigner(TString option = ""){
        m_w->Draw(option);
    }
    //return the the normalization constant
    double const NormalizationConstant(){
        return m_const;
    }

    //W x W
    double const WignerForItself(){
        TString wtimesw_name = m_name + "timesItslef";
        TF2 *wTimesw = new TF2(wtimesw_name,wignerSource2,m_minX,m_maxX,m_minP,m_maxP,5);
        wTimesw->SetParameter(0, m_const);
        wTimesw->SetParameter(1, m_r0);
        wTimesw->SetParameter(2, m_k);
        wTimesw->SetParameter(3, m_hcut);
        wTimesw->SetParameter(4, m_ndim);


        return myintegral(wTimesw)*TMath::Power(4*m_hcut*TMath::Pi(),m_ndim);
    }

    void SetR(float r0){
        m_r0 = r0;
        m_w->SetParameter(1, m_r0); 
        m_const = normalization();
        m_w->SetParameter(0, m_const);
    }
    void SetK(float k){
        m_k = k;
        m_w->SetParameter(2, m_k); 
        m_const = normalization();
        m_w->SetParameter(0, m_const);
    }
    void SetRK(float r0, float k){
        m_r0 = r0;
        m_k = k;
        m_w->SetParameter(1, m_r0); 
        m_w->SetParameter(0, m_const);
        m_const = normalization();
        m_w->SetParameter(0, m_const);

    }

};



void WignerSource(){
    auto fw = wignerfunction(5., 50., 3, 0.,50.,0,300, "w");
    std::cout<<"Normalization: "<<fw.NormalizationConstant()<<"\n";
    std::cout<<"W x W = "<<fw.WignerForItself()<<"\n";
    fw.SetRK(5., 200);
    std::cout<<"Normalization: "<<fw.NormalizationConstant()<<"\n";
    std::cout<<"W x W = "<<fw.WignerForItself()<<"\n";
    fw.DrawWigner("SURF2");
}