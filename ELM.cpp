#include "ELM.h"

CEl::CEl(const int& NEl, const int& DegreeApproxPol, const int& NMa, tvMa& vMa, const int& NNo, const V1I& vNNo, tvNo& vNo, const int& NQuad, tvQuadratura& Qua, CCo& Co)
{
  this->NEl=NEl;
  this->dimension=2;
  this->DegreeApproxPol=DegreeApproxPol;
  this->NNo=NNo;
  dDegreeApproxPol=DegreeApproxPol-1;
  NNEl1m=c_NNEl1m();
  for (int a=0; a<this->NNo; a++)
  {
	No.push_back(&vNo[vNNo[a]]);
	no.push_back( vNo[vNNo[a]].r_NEl());
	No[a]->w_FP(1);
  }
  this->NMa=NMa;
  this->Ma=&vMa[NMa];
  this->NQuad=NQuad;
  this->Quad=&Quad[NQuad];
  this->Co=&Co;
  Xsi.resize(NNo);
  FF.resize(NNo);
  no_E.resize(NNo);
  no_S.resize(NNo);
  no_ue.resize(NNo,0.0);
  dFF.resize(NNEl1m);
  for (int i = 0; i < NNo; ++i)
  {
	Xsi[i].resize(dimension,0.0);
	FF[i].resize(NNo,0.0);
	no_E[i].resize(dimension*dimension,0.0);
	no_S[i].resize(dimension*dimension,0.0);
  }
  for (int i = 0; i < NNEl1m; ++i)
  {
	dFF[i].resize(NNo,V1D(dimension,0.0));
  }
  Fint.resize(dimension*NNo,0.0);
  Hessiana.resize(dimension*NNo,V1D(dimension*NNo,0.0));
//  Mass.resize(dimension*NNo,V1D(dimension*NNo,0.0));
  
  PHI.resize(r_Quad()->r_NumberQuad());
  DPHI.resize(r_Quad()->r_NumberQuad());
  X.resize(r_Quad()->r_NumberQuad());
  A0.resize(r_Quad()->r_NumberQuad());
  J0.resize(r_Quad()->r_NumberQuad(),0.0);
  invA0.resize(r_Quad()->r_NumberQuad());
  Y.resize(r_Quad()->r_NumberQuad());
  A1.resize(r_Quad()->r_NumberQuad());
  J1.resize(r_Quad()->r_NumberQuad(),0.0);
  A.resize(r_Quad()->r_NumberQuad());
  C.resize(r_Quad()->r_NumberQuad());
  E.resize(r_Quad()->r_NumberQuad());
  S.resize(r_Quad()->r_NumberQuad());
  CC.resize(r_Quad()->r_NumberQuad());
  ue.resize(r_Quad()->r_NumberQuad(),0.0);
  fint.resize(r_Quad()->r_NumberQuad());
  h.resize(r_Quad()->r_NumberQuad());
  //m.resize(r_Quad()->r_NumberQuad());
  for (int i = 0; i < r_Quad()->r_NumberQuad(); ++i)
  {
	PHI[i].resize(NNo,0.0);
	DPHI[i].resize(NNo,V1D(dimension,0.0));
	X[i].resize(dimension,0.0);
	A0[i].resize(dimension,V1D(dimension,0.0));
	invA0[i].resize(dimension,V1D(dimension,0.0));
	Y[i].resize(dimension,0.0);
	A1[i].resize(dimension,V1D(dimension,0.0));
	A[i].resize(dimension,V1D(dimension,0.0));
	C[i].resize(dimension,V1D(dimension,0.0));
	E[i].resize(dimension*dimension,0.0);
	S[i].resize(dimension*dimension,0.0);
	CC[i].resize(dimension*dimension,V1D(dimension*dimension,0.0));
	fint[i].resize(dimension*NNo,0.0);
	h[i].resize(dimension*NNo,V1D (dimension*NNo,0.0));
	//m[i].resize(dimension*NNo,V1D (dimension*NNo,0.0));
  }
  Xsi=c_Xsi(NNo,DegreeApproxPol,dimension);
  FF =c_FF (NNo,DegreeApproxPol,Xsi);
  dFF=c_dFF(NNo,DegreeApproxPol,dimension,dDegreeApproxPol,NNEl1m,FF);
  Ue=0.0;
  for (int i=0; i<r_Quad()->r_NumberQuad(); i++)
  {
	PHI[i] =calc_PHI (dDegreeApproxPol,NNo,dimension,r_Quad()->r_Xsi(i),FF);
	DPHI[i]=calc_DPHI(dDegreeApproxPol,NNo,dimension,r_Quad()->r_Xsi(i),dFF);
    A0[i]=calc_An(DPHI[i],No,0);
    J0[i]=Det(A0[i]);
	X[i]=calc_Xn(dimension,PHI[i],No,0);
	invA0[i]=inverse(A0[i]);  //www
  }
}




void CEl::ite()
{
  double beta_=1.0;
  double deltat_=1.0;
    
    
    
    
  for (int i=0; i<dimension*NNo; i++)
  {
    Fint[i]=0.0;
	for (int j=0; j<dimension*NNo; j++)
    {
	  Hessiana[i][j]=0.0;
	  //Mass[i][j]=0.0;
    }
  }
  Ue=0.0;
  for (int i=0; i<r_Quad()->r_NumberQuad(); i++)
  {
	Y[i]=calc_Xn(dimension,  PHI[i],No,1);
    A1[i]=calc_An(DPHI[i],No,1); // Gradiente de tranformação final
    J1[i]=Det(A1[i]);                    
	A[i]=MM(A1[i],invA0[i]);
    C[i]=MtM(A[i],A[i]);    
    E[i]=calc_Ei(C[i]); //executa a cada iteracao para cada ponto de gauss e armazena Si=due_dE, e CCi=d2ue_dEdE
	CC[i]=calc_CC(dimension,Ma); // Calculo do tensor contitutivo de Piola-Kirchhoff de segunda espécie;
	S[i]=calc_S(dimension,E[i],CC[i]); // Calculo da Tensões de Green;
	calc_Fint_H_M(PHI[i],DPHI[i],invA0[i],A[i],S[i],CC[i],fint[i],h[i],m[i]);   // Calculo do vetor de forças internas e da matriz Hesiana
      
      
      
      
      
      
      
    ue[i]=c_ue(E[i],Ma);
	Ue=Ue+ue[i]*J0[i]*r_Quad()->r_W(i);
	for (int a=0; a<NNo; a++)
    {
	  for (int z=0; z<dimension; z++)
      {
		Fint[dimension*a+z]+=fint[i][dimension*a+z]*J0[i]*r_Quad()->r_W(i);
		for (int l=0; l<NNo; l++)
        {
		  for (int k=0; k<dimension; k++)
          {
			Hessiana[dimension*a+z][dimension*l+k]+=h[i][dimension*a+z][dimension*l+k]*J0[i]*r_Quad()->r_W(i);
			if (z==k) Mass[d*a+z][dimension*l+k] += (m[i][dimension*a+z][dimension*l+k]/(beta_ * deltat_ * deltat_))*J0[i]*r_Quad()->r_W(i);/**/ //mass contribution
          }
        }
      }
	}
  }
  for (int l=0; l<NNo; l++)
  {
	cpto(dDegreeApproxPol,NNo,dimension,Xsi[l],FF,dFF,No,Ma,no_E[l],no_S[l],no_ue[l]);
  }
  Imp_Rel_Geral();
}





void cpto(const int& dDegreeApproxPol,const int& NNo,const int& dimension,const V1D XsiNo,const V2D& FF,const V3D& dFF,  std::vector< CNo * > No, CMa *Ma, V1D& E, V1D& S, double& ue)
{
  V1D    PHI;  ///< Valor das funções de forma nos pontos de integração
  V2D    DPHI; ///< Valor da primeira derivada das funcoes de forma nos pontos de integração
  V1D    X;    ///< Posição do ponto de integração na configuração inicial
  V2D    A0;   ///< Gradiente de tranformação inicial
  double J0;   ///< Determinante do Gradiente de tranformação inicial
  V2D    invA0;    ///< Matriz inversa da Gradiente de tranformação inicial
  V1D    Y;    ///< Posição do ponto de integração na configuração final
  V2D    A1;   ///< Gradiente de tranformação final
  double J1;   ///< Determinante do Gradiente de tranformação final
  V2D    A;    ///< Gradiente de transformação total
  V2D    C;    ///< Tensor de alongamento ou estiramento a direita de Cauchy-Green
  V2D    CC;   ///< Tensor constitutivo de Piola Kirchhoff de segunda espécie
  V1D    fint; ///< Vetor de forças específica* internas
  V2D    h;    ///< Matrix Hesiana específica*
  V2D    m;    ///< Matrix de massa específica*

  PHI.resize(NNo,0.0);
  DPHI.resize(NNo,V1D(dimension,0.0));
  X.resize(dimension,0.0);
  A0.resize(dimension,V1D(dimension,0.0));
  J0=0.0;
  invA0.resize(dimension,V1D(dimension,0.0));
  Y.resize(dimension,0.0);
  A1.resize(dimension,V1D(dimension,0.0));
  J1=0.0;
  A.resize(dimension,V1D(dimension,0.0));
  C.resize(dimension,V1D(dimension,0.0));
  CC.resize(dimension*dimension,V1D(dimension*dimension,0.0));
  fint.resize(dimension*NNo,0.0);
  h.resize(dimension*NNo,V1D (dimension*NNo,0.0));
  m.resize(dimension*NNo,V1D (dimension*NNo,0.0));
//  std::cout << "Xsi " << XsiNo[0] << " " << XsiNo[1] << " " << XsiNo[2] << std::endl;
//  std::cout << "FF"<< std::endl; for (int jj=0; jj<NN; jj++) {   for (int ii=0; ii<NN; ii++) { std::cout << FF[jj][ii] << " "; } std::cout<<std::endl; }
  PHI =calc_PHI (dDegreeApproxPol,NNo,dimension,XsiNo,FF);//std::cout << "PHI"<< std::endl; for (int ii=0; ii<PHI.size(); ii++) { std::cout << PHI[0] << " "; } std::cout<<std::endl;
  DPHI=calc_DPHI(dDegreeApproxPol,NNo,dimension,XsiNo,dFF);//std::cout << "DPHI"<< std::endl; for (int jj=0; jj<d; jj++) {   for (int ii=0; ii<DPHI.size(); ii++) { std::cout << DPHI[ii][jj] << " "; } std::cout<<std::endl; }
  A0=calc_An(DPHI,No,0);
  J0=Det(A0);
  X=calc_Xn(dimension,PHI,No,0);//  std::cout << "X " << X[0] << " " << X[1] << " " << X[2] << std::endl;
  invA0=inverse(A0);
  Y=calc_Xn(dimension,  PHI,No,1); //  std::cout << "Y " << Y[0] << " " << Y[1] << " " << Y[2] << std::endl;
  A1=calc_An(DPHI,No,1); // Gradiente de tranformação final
  J1=Det(A1);                    
  A=MM(A1,invA0); // std::cout << "A " << A[0][0] << " " << A[0][1] << " " << A[0][2] << " " << A[1][0] << " " << A[1][1] << " " << A[1][2] << " " << A[2][0] << " " << A[2][1] << " " << A[2][2] << std::endl;
  C=MtM(A,A);  //std::cout << "C " << C[0][0] << " " << C[0][1] << " " << C[0][2] << " " << C[1][0] << " " << C[1][1] << " " << C[1][2] << " " << C[2][0] << " " << C[2][1] << " " << C[2][2] << std::endl; 
  E=calc_Ei(C); // std::cout << "E " << E[0][0] << " " << E[0][1] << " " << E[0][2] << " " << E[1][0] << " " << E[1][1] << " " << E[1][2] << " " << E[2][0] << " " << E[2][1] << " " << E[2][2] << std::endl; 
  CC=calc_CC(dimension,Ma); // Calculo do tensor contitutivo de Piola-Kirchhoff de segunda espécie;
  S=calc_S(dimension,E,CC); // Calculo da Tensões de Green;
  calc_Fint_H_M(PHI,DPHI,invA0,A,S,CC,fint,h,m);   // Calculo do vetor de forças internas e da matriz Hesiana
  ue=c_ue(E,Ma);
}


CEl::~CEl()
{
  No.clear();
  no.clear();
  Xsi.clear();
  FF.clear();
  dFF.clear();
  Fint.clear();
  Hessiana.clear();
  PHI.clear();
  DPHI.clear();
  X.clear();
  A0.clear();
  J0.clear();
  invA0.clear();
  Y.clear();
  A1.clear();
  J1.clear();
  A.clear();
  C.clear();
  E.clear();
  CC.clear();
  S.clear();
  ue.clear();
  fint.clear();
  h.clear();
  //m.clear();
}


int CEl::c_NNEl1m()
{
  int nn1m;
  nn1m=(dDegreeApproxPol+1)*((dDegreeApproxPol+1)+1)/2; // numero de nos do elemento triangular com grau do polinomio uma vez menor
  return nn1m;
}

V2D c_Xsi(const int& NNo, const int& DegreeApproxPol, const int& dimension)
{
  V2D Xsi; Xsi.resize(NNo, V1D (dimension,0.0));
// NN que é dado pode também ser calculado por: NN=(Pe+1)*(2+Pe)/2;		
// Calcula as Coordenadas Adimensionais dos Nós
  for (int I=0; I<(DegreeApproxPol+1); I++)
  {
	for (int J=0; J<(I+1); J++)
    {
      int K;
      K=I+((I-1)*I)/2+J;
	  Xsi[K][0]=(I/(1.0*DegreeApproxPol))-(J)/(1.0*DegreeApproxPol);				// coordenadas adimensionais dos nos
	  Xsi[K][1]=(J)/(1.0*DegreeApproxPol);
    }
  }
  return Xsi;
}



V2D c_FF(const int& NNo, const int& DegreeApproxPol, const V2D& Xsi)
{
  V2D MAT;        MAT.resize(NNo,V1D(NNo,0.0));
  V2D FF;         FF.resize(NNo,V1D(NNo,0.0));
// NN que é dado pode também ser calculado por: NN=(Pe+1)*(2+Pe)/2;		
// Calcula a Matriz a ser Invertida para a determinacao dos Coeficientes das Funcoes de forma
  for (int L=0; L<NNo; L++)
  {
	for (int I=0; I<(DegreeApproxPol+1); I++)
    {
      for (int J=0; J<(I+1); J++) 
	  {
        int K;
        K=I+((I-1)*I)/2+J;
        MAT[L][K]=(pow(Xsi[L][0],(I-J)))*(pow(Xsi[L][1],(J)));		// matriz a ser invertida para calculo dos coeficientes...
      }
    }
  }
// Calcula a Matriz dos Coeficientes das Funcoes de forma
  FF=inverse(MAT);         // ANTERIORMENTE ENUMERADO NA FUNÇÃO LAPAC
  for (int I=0; I<NNo; I++)
  {
	for (int J=0; J<NNo; J++)
    {
	  if (fabs(FF[I][J])<0.0000000001) FF[I][J]=0.0;
    }
  }
  return FF;
}

V3D c_dFF(const int& NNo, const int& DegreeApproxPol, const int& dimension, const int& dDegreeApproxPol, const int& NNEl1m, const V2D& FF)
{
  V3D dFF;           dFF.resize(NNEl1m,V2D(NNo,V1D(dimension,0.0)));
// Calcula a Matriz dos Coeficientes das Derivadas direcionais das Funcoes de forma
  for (int L=0; L<NNo; L++)
  {
	for (int I=0; I<(dDegreeApproxPol); I++)
    {
      for (int J=0; J<(I+1); J++) 
      {
		int K;
        K=I+((I-1)*I)/2+J;
		dFF[K][L][0]=FF[K+I+1][L]*(I+1-(J));			// coeficientes das derivadas das funÁıes de forma
		dFF[K][L][1]=FF[K+I+2][L]*(J+1);
      }
    }
  }
  return dFF;
}


V1D calc_PHI(const int& DegreeApproxPol, const int& NNo, const int dimension, const V1D& Xsi, const V2D& FF)
{
  V1D PHIi;
  PHIi.resize(NNo,0.0);
 // calcula os Valores das funcoes de forma no ponto de Integracao atual(PHI)
  for (int L=0; L<NNo; L++)
  {
	PHIi[L]=0.0;
	for (int I=0; I<(DegreeApproxPol+1); I++)
    {
      for (int J=0; J<(I+1); J++) 
      {
        int K;    
		K=I+((I-1)*I)/2+J;
		PHIi[L]+=FF[K][L]*(pow(Xsi[0],(I-J)))*(pow(Xsi[1],(J)));	// valores das funcoes de forma no ponto com...
      }															// ...coordenadas (COORD1,COORD2)
    }
    if (fabs(PHIi[L])<0.0000000001) PHIi[L]=0.0;
  }
  return PHIi;
}


V2D calc_DPHI(const int& DegreeApproxPol, const int& NNo, const int dimension, const V1D& Xsi, const V3D& dFF)
{
  V2D DPHIi;
  DPHIi.resize(NNo,V1D(dimension,0.0));
  for (int M=0; M<2; M++)
  {
	for (int L=0; L<NNo; L++)
    {
	  for (int I=0; I<(DegreeApproxPol); I++)
      {
        for (int J=0; J<(I+1); J++) 
        {
          int K;    
          K=I+((I-1)*I)/2+J;
		  DPHIi[L][M]+=dFF[K][L][M]*(pow(Xsi[0],(I-J)))*(pow(Xsi[1],(J)));		// valores das derivadas das funÁıes de forma no...//...ponto (COORD1,COORD2)
        }
      }
      if (fabs(DPHIi[L][M])<0.0000000001) DPHIi[L][M]=0.0;
    }
  }
  return DPHIi;
}


V1D calc_Xn(const int& dimension, const V1D& PHI, std::vector< CNo * > No, const int estado)
{
  V1D xi; xi.resize(dimension,0.0);
  int NN=PHI.size();
  for (int m=0; m<dimension; m++)
  {
	xi[m]=0.0;
	for (int L=0; L<NN; L++)
    {
      xi[m]+=PHI[L]*No[L]->r_X(estado,m);
    }
  }
  return xi;
}

V2D calc_An(const V2D& DPHI, std::vector< CNo * > No, const int estado)
{
  int NN=DPHI.size();
  int d=DPHI[0].size();
  V2D Ani; Ani.resize(d, V1D (d,0.0));
  for (int m=0; m<d; m++) 
  {
    for (int n=0; n<d; n++) 
    {
      Ani[m][n]=0.0;
      for (int L=0; L<DPHI.size(); L++)
      {
        Ani[m][n]+=DPHI[L][n]*No[L]->r_X(estado,m);
      }
    }
  }
  return Ani;
}


V1D calc_Ei(const V2D& Ci)
{
  V1D Ei; Ei.resize(Ci.size()*Ci[0].size(),0.0); 
  Ei[0]=(Ci[0][0]-1.0)/2.0;
  Ei[1]=(Ci[1][1]-1.0)/2.0;
  Ei[2]=(Ci[0][1]    )/2.0;
  Ei[3]=(Ci[1][0]    )/2.0;
  return Ei;
}


V2D calc_CC(const int& dimension, CMa *Ma)
{
  V2D CCi;
  CCi.resize(dimension*dimension,V1D(dimension*dimension,0.0));
  CCi[0][0]=(2*Ma->r_LameM()+Ma->r_LameL());
  CCi[0][1]=   Ma->r_LameL()               ;
  CCi[1][0]=   Ma->r_LameL()               ;
  CCi[1][1]=(2*Ma->r_LameM()+Ma->r_LameL());
//verificar direito estes termos:
  CCi[2][2]=Ma->r_LameM();
  CCi[2][3]=Ma->r_LameM();
  CCi[3][2]=Ma->r_LameM();
  CCi[3][3]=Ma->r_LameM();
  return CCi;
}


V1D calc_S(const int& dimension, const V1D& Ei, V2D& CCi)
{
  V1D Si;
  Si.resize(dimension*dimension,0.0);
  Si[0]=CCi[0][0]*Ei[0]+CCi[0][1]*Ei[1];
  Si[1]=CCi[1][0]*Ei[0]+CCi[1][1]*Ei[1];
  Si[2]=CCi[2][2]*Ei[2]+CCi[2][3]*Ei[3];
  Si[3]=CCi[3][2]*Ei[2]+CCi[3][3]*Ei[3];
  return Si;
}


double c_ue(const V1D& E, CMa *Ma)
{
  double ue;
  ue=(1.0/2.0)*
      (
	   (2.0*Ma->r_LameM()+Ma->r_LameL())*(E[0]*E[0]+E[1]*E[1]    )+
	   (                  Ma->r_LameL())*(E[0]*E[1]              )+
	   (    Ma->r_LameM()              )*((E[2]+E[3])*(E[2]+E[3]))
	  );
  return ue;
}


void calc_Fint_H_M(const V1D& PHI, const V2D& DPHI, const V2D& Di, const V2D& Ai, const V1D& Si, const V2D& CCi, V1D& finti, V2D& hi, V2D& mi)
{
  int NN=DPHI.size();
  int d=DPHI[0].size();
  V4D dA1_dY;          dA1_dY.resize(NN,V3D(d,V2D(d,V1D(d,0.0))));
  V4D dA_dY;           dA_dY.resize (NN,V3D(d,V2D(d,V1D(d,0.0))));
  V1D dE_dYaz;         dE_dYaz.resize(d*d,0.0); 
  V1D dE_dYlk;         dE_dYlk.resize(d*d,0.0);
  V1D d2E_dYdY;        d2E_dYdY.resize(d*d,0.0);
  V1D CCi_dE_dYlk;     CCi_dE_dYlk.resize(d,0.0);
  V2D AdA_dY;          AdA_dY.resize(d,V1D(d,0.0));
  V2D dA_dY_dA_dY_1;   dA_dY_dA_dY_1.resize(d,V1D(d,0.0)); 
  V2D dA_dY_dA_dY_2;   dA_dY_dA_dY_2.resize(d,V1D(d,0.0));
    
    
    double density=1.0;
    
    
    
  for (int L=0; L<NN; L++)
  {
	for (int m=0; m<d; m++)
	{
	  for (int n=0; n<d; n++)
	  {
		dA1_dY[L][m][m][n]=DPHI[L][n];
	  }
	}
  }
  for (int a=0; a<NN; a++)
  {
    for (int z=0; z<d; z++)
    {
	  dA_dY[a][z]=MM(dA1_dY[a][z],Di);
	}
  }

  //Calcula Matriz hessiana e fint no ponto de integração
  for (int a=0; a<NN; a++)
  {
	for (int z=0; z<d; z++)
	{
	  AdA_dY=MtM(Ai,dA_dY[a][z]);
	  dE_dYaz[0]=(1.0/2.0)*(AdA_dY[0][0]+AdA_dY[0][0]);
	  dE_dYaz[1]=(1.0/2.0)*(AdA_dY[1][1]+AdA_dY[1][1]);
	  dE_dYaz[2]=(1.0/2.0)*(AdA_dY[1][0]+AdA_dY[0][1]);
	  dE_dYaz[3]=(1.0/2.0)*(AdA_dY[0][1]+AdA_dY[1][0]);
	  finti[2*a+z]=Si[0]*dE_dYaz[0]+Si[1]*dE_dYaz[1]+Si[2]*dE_dYaz[2]+Si[3]*dE_dYaz[3];
      for (int l=0; l<NN; l++)
      {
        for (int k=0; k<d; k++)
        {
          AdA_dY=MtM(Ai,dA_dY[l][k]);
          dE_dYlk[0]=(1.0/2.0)*(AdA_dY[0][0]+AdA_dY[0][0]);
          dE_dYlk[1]=(1.0/2.0)*(AdA_dY[1][1]+AdA_dY[1][1]);
          dE_dYlk[2]=(1.0/2.0)*(AdA_dY[1][0]+AdA_dY[0][1]);
          dE_dYlk[3]=(1.0/2.0)*(AdA_dY[0][1]+AdA_dY[1][0]);
          dA_dY_dA_dY_1=MtM(dA_dY[a][z],dA_dY[l][k]);
          dA_dY_dA_dY_2=MtM(dA_dY[l][k],dA_dY[a][z]);
          d2E_dYdY[0]=(1.0/2.0)*(dA_dY_dA_dY_1[0][0]+dA_dY_dA_dY_2[0][0]);
          d2E_dYdY[1]=(1.0/2.0)*(dA_dY_dA_dY_1[1][1]+dA_dY_dA_dY_2[1][1]);
          d2E_dYdY[2]=(1.0/2.0)*(dA_dY_dA_dY_1[0][1]+dA_dY_dA_dY_2[0][1]);
          d2E_dYdY[3]=(1.0/2.0)*(dA_dY_dA_dY_1[1][0]+dA_dY_dA_dY_2[1][0]);
          CCi_dE_dYlk=MV(CCi,dE_dYlk);
          hi[2*a+z][2*l+k]=dE_dYaz[0]*CCi_dE_dYlk[0]+
                           dE_dYaz[1]*CCi_dE_dYlk[1]+
                           dE_dYaz[2]*CCi_dE_dYlk[2]+
                           dE_dYaz[3]*CCi_dE_dYlk[3]+
                           (Si[0]*d2E_dYdY[0]+
                            Si[1]*d2E_dYdY[1]+
                            Si[2]*d2E_dYdY[2]+
                            Si[3]*d2E_dYdY[3]);
            
            if (z==k) mi[2*a+z][2*l+k] = density * PHI[a] * PHI[l]; //mass contribution
            
            
            
            
        }
      }
    }
  }
    
    
    
    
    
    
    
    
    
    
    
    
    
}


void Le_Elementos(tvEl& El, tvNo& No, tvMa& Ma, tvQuadratura& Qua, CCo& Co)
{
  int i,nEl;
  char s[1000];
  std::ifstream Ent; // para leitura
  std::string NAr=Co.r_ArqELM();
  Ent.open(NAr.c_str());
  Ent.getline(s,1000);
  Ent >>  nEl; Ent.getline(s,1000);
  Ent.getline(s,1000);
  for (i=0; i<nEl; i++) 
  {
    V1I vNNo;
    int nno;  ///< Numero de nos do elemento;
    int D=2;    ///< Dimensão do elemento: 1D, 2D ou 3D;
    int Pe;   ///< Aproximacao fisica do elemento: 1=linear 2=quadratico ...;
    int ma;   ///< Indice do material que compoe o elemento;
    int qua;  ///< Indice da quadratura a ser usada para integrar o elemento;
    int I;    ///< Indice do elemento;
    Ent >> I >> Pe >> ma >> qua >> nno;
    vNNo.resize(nno,0);
    for (int j=0; j<nno; j++) 
    {
      Ent >> vNNo[j];
    }
    Ent.getline(s,1000);
    CEl LeEl(i, Pe, ma ,Ma, nno, vNNo, No, qua, Qua, Co);
    El.push_back(LeEl);
  }
  Ent.close();
}

void CEl::Imp_Rel_Geral()
{
  std::string NArq1; NArq1.append(Co->r_ArqREL());
  std::fstream Sai1;
  Sai1.open(NArq1.c_str(),std::fstream::in | std::fstream::out | std::fstream::app);     ///< Dados Calculados
  Sai1 << "Iteração =" << "X" << std::endl;
  Sai1 << "Elemento =" << NEl << std::endl;
  Sai1 << "GRAU DO POLINOMIO =" << DegreeApproxPol << std::endl;
  Sai1 << "NUMERO DE NOS =" << NNo << std::endl;
  Sai1 << "Coordenadas Adimensionais dos Nós" << std::endl;
  for (int I=0; I<Xsi.size(); I++) { Sai1 << I << " "; for (int J=0; J<Xsi[0].size(); J++) { Sai1 << "Xsi[" << I << "][" << J << "]=" << Xsi[I][J] << "   "; } Sai1 << std::endl; }
  Sai1 << "Matriz dos Coeficientes das Funcoes de forma" << std::endl;
  for (int I=0; I<NNo; I++) { for (int J=0; J<NNo; J++) { Sai1 << FF[I][J] << " "; } Sai1 << std::endl; }
  Sai1 << "Matriz dos Coeficientes das Derivadas direcionais das Funcoes de forma" << std::endl;
  Sai1 << "C'0" << std::endl; for (int I=0; I<NNEl1m; I++) { for (int J=0; J<NNo; J++) { Sai1 << dFF[I][J][0] << " "; } Sai1 << std::endl; }
  Sai1 << "C'1" << std::endl; for (int I=0; I<NNEl1m; I++) { for (int J=0; J<NNo; J++) { Sai1 << dFF[I][J][1] << " "; } Sai1 << std::endl; }
  Sai1 << std::endl;
  for (int i=0; i<r_Quad()->r_NumberQuad(); i++)
  {
    Sai1 << "Ponto de integracao No. =" << i << std::endl;
	Sai1 << "Valores das funcoes de forma no ponto de Integracao atual(PHI) :" << std::endl;
	for (int L=0; L<FF[0].size(); L++) { Sai1 << PHI[i][L] << " "; } Sai1 << std::endl;
    Sai1 << "Valores das derivadas direcionais das funcoes de forma no ponto de Integracao atual(DPHI):" << std::endl;
	for (int M=0; M<dimension; M++) { for (int L=0; L<dFF[0].size(); L++) { Sai1 << DPHI[i][L][M] << " "; } Sai1 << std::endl; }
	Sai1 << " X= " << std::endl; for (int m=0; m<dimension; m++) { Sai1 << X[i][m] << " "; } Sai1 << std::endl;
    Sai1 << " J0= " << J0[i] << std::endl;
	Sai1 << " Y= " << std::endl; for (int m=0; m<dimension; m++) { Sai1 << Y[i][m] << " "; } Sai1 << std::endl;
	Sai1 << "Matriz A1= " << std::endl; for (int m=0; m<dimension; m++) { for (int n=0; n<dimension; n++) { Sai1 << A1[i][m][n] << " "; } Sai1 << std::endl; }
    Sai1 << " J1= " << J1[i] << std::endl;
	Sai1 << " A= " << std::endl; for (int m=0; m<dimension; m++) { for (int n=0; n<dimension; n++) { Sai1 << A[i][m][n] << " "; } Sai1 << std::endl; }
	Sai1 << " C= " << std::endl; for (int m=0; m<dimension; m++) { for (int n=0; n<dimension; n++) { Sai1 << C[i][m][n] << " "; } Sai1 << std::endl; }
	Sai1 << " E= " << std::endl; for (int m=0; m<dimension*dimension; m++) { Sai1 << E[i][m] << " "; } Sai1 << std::endl;
	Sai1 << " S= " << std::endl; for (int m=0; m<dimension*dimension; m++) { Sai1 << S[i][m] << " "; } Sai1 << std::endl;
	Sai1 << " CC= " << std::endl; for (int m=0; m<dimension*dimension; m++) { for (int n=0; n<dimension*dimension; n++) { Sai1 << CC[i][m][n] << " "; } Sai1 << std::endl; }
	Sai1 << " fint= "; for (int j=0; j<dimension*NNo; j++) { Sai1 << fint[i][j] << " "; } Sai1 << std::endl;
	Sai1 << " h= "<< std::endl; for (int a=0; a<dimension*NNo; a++) { for (int j=0; j<dimension*NNo; j++) { Sai1 << h[i][a][j] << " "; } Sai1 << std::endl; }
	Sai1 << " ue[" << i << "]=" << ue[i] << " J0=" << J0[i] << " W[i]=" << r_Quad()->r_W(i) <<  std::endl;
    Sai1 << std::endl;
  }
  Sai1 << " Ue= " << Ue <<  std::endl;
  Sai1 << " Fint= " << std::endl; for (int i=0; i<dimension*NNo; i++) { Sai1 << Fint[i] << " "; } Sai1 << std::endl;
  Sai1 << " H= " << std::endl; for (int i=0; i<dimension*NNo; i++) { for (int j=0; j<dimension*NNo; j++) { Sai1 << Hessiana[i][j] << " "; } Sai1 << std::endl; }
  Sai1 << " M= " << std::endl; for (int i=0; i<dimension*NNo; i++) { for (int j=0; j<dimension*NNo; j++) { Sai1 << Mass[i][j] << " "; } Sai1 << std::endl; }
  Sai1 << std::endl;
  Sai1 << std::endl;
  Sai1.close();
}
