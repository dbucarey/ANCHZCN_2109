//#########################################################
GLOBALS_SECTION
//#########################################################

 #include <admodel.h>
 #include <stdio.h>
 #include <time.h>
 #include <limits>
 time_t start,finish;
 long hour,minute,second;
 double elapsed_time;

 double eps=std::numeric_limits<double>::epsilon();

 #include <string.h>
 #undef veo
 #undef depuro
 #define veo(object) cout << #object "\n" << object << endl;
 #define depuro(object) cout << #object "\n" << object << endl; exit(1);


//#########################################################
TOP_OF_MAIN_SECTION
//#########################################################

 time(&start);
 arrmblsize = 90000000;
 gradient_structure::set_GRADSTACK_BUFFER_SIZE(1.e7);
 gradient_structure::set_CMPDIF_BUFFER_SIZE(1.e7);
 gradient_structure::set_MAX_NVAR_OFFSET(5000);

//#########################################################
DATA_SECTION
//#########################################################
 init_int ntime
 init_int nedades
 init_int ntallas
 init_vector edades(1,nedades)
 init_vector Tallas(1,ntallas)
 init_matrix data(1,ntime,1,13)///Yrs, desem, cv1, CPUEI, cv2, CPUEA, cv3, Bcru, cv4, MPH, cv5, nmf	nmc
 init_matrix Ctot(1,ntime,1,ntallas)
 init_matrix Ncru(1,ntime,1,ntallas)
 init_vector msex(1,ntallas)
 init_vector Wmed(1,ntallas)

//========================================================
//LEER controles y opciones
//!! ad_comm::change_datafile_name(".ctl");
//========================================================

//  Coeficiente de variacion de los desvios Rt, No
 init_vector cvar(1,2) //cv_Rt y cv_No
 init_vector dt(1,3)//dt_desove, cpue, reclan
//-------------------------------------------------------------
// parametros de crecimiento,alfa y beta relacion stock recluta,
// mortalidad natural,y stepness
//-------------------------------------------------------------
 init_vector Par_bio(1,7)//Loo,k,Lo,alfa,beta,M,h
 init_vector cv_priors(1,7)//cv prior si Par_bio son estimados
 init_vector fases_bio(1,7)//fases de estimacion Par_bio

//-------------------------------------
// calculos preliminares
//-------------------------------------
  //parametros
  number log_Linfprior
  number log_kprior
  number log_Loprior
  number log_alfa_prior
  number log_beta_prior
  number log_M_prior
  number log_h_prior

  !! log_Linfprior  = log(Par_bio(1));
  !! log_kprior     = log(Par_bio(2));
  !! log_Loprior    = log(Par_bio(3));
  !! log_alfa_prior = log(Par_bio(4));
  !! log_beta_prior = log(Par_bio(5));
  !! log_M_prior    = log(Par_bio(6));
  !! log_h_prior    = log(Par_bio(7));
//-------------------------------------
  //fases de estimacion parametros
  number opt_Linf
  number opt_k
  number opt_Lo
  number opt_alfa
  number opt_beta
  number opt_M
  number opt_h

  !! opt_Linf = fases_bio(1);
  !! opt_k    = fases_bio(2);
  !! opt_Lo   = fases_bio(3);
  !! opt_alfa = fases_bio(4);
  !! opt_beta = fases_bio(5);
  !! opt_M    = fases_bio(6);
  !! opt_h    = fases_bio(7);

//--------------------------------------
// Capturabilidad y cv capturabilidad
//--------------------------------------
 init_number prior_qf//q flota industrial
 init_number prior_qfa//q flota artesanal
 init_number prior_qc//q reclan
 init_number prior_qmph//q mpdh
 init_number cv_qf//cv q cpue industrial
 init_number cv_qfa//cv q cpue artesanal
 init_number cv_qcru//cv q reclan
 init_number cv_qmph//cv q mph


  //calculo preliminar
  number log_prior_qf
  number log_prior_qfa
  number log_prior_qc
  number log_prior_qmph
  !! log_prior_qf = log(prior_qf);
  !! log_prior_qfa = log(prior_qfa);
  !! log_prior_qc = log(prior_qc);
  !! log_prior_qmph = log(prior_qmph);


//--------------------------------------
//SELECTIVIDADES - A50 y rango
//--------------------------------------
 init_vector parf(1,2) //valores de partida A50 (edad) y rango
 init_vector parc(1,2)

 number log_L50priorf
 number log_spriorf
 number log_L50priorc
 number log_spriorc
 !! log_L50priorf = log(parf(1));
 !! log_spriorf   = log(parf(2));
 !! log_L50priorc = log(parc(1));
 !! log_spriorc   = log(parc(2));

//--------------------------------------
// BLOQUES DE SELECTIVIDAD
//--------------------------------------
 init_int    nbloques1 //numero de bloques de selectividad flota
 init_vector ybloques1(1,nbloques1)//año de inicio bloque de selectividad flota
 init_int    nbloques2//numero de bloques de selectividad RECLAN
 init_vector ybloques2(1,nbloques2)//año de inicio bloque de selectividad RECLAN

//--------------------------------------
// BLOQUES DE CAPTURABULIDAD FFLOTA y CRUCEROS
//--------------------------------------
 init_int    nqbloques//numero de bloques de Capturabilidad cpue industrial
 init_vector yqbloques(1,nqbloques)//año de inicio bloque de Capturabilidad cpue industrial

 init_int    nqabloques//bloques//numero de bloques de Capturabilidad cpue artesanal
 init_vector yqabloques(1,nqabloques)//año de inicio bloque de Capturabilidad cpue artesanal

 init_int    nqbloquesc//numero de bloques de Capturabilidad RECLAN
 init_vector yqbloquesc(1,nqbloquesc)//año de inicio bloque de Capturabilidad RECLAN

 init_int    nqbloquesmph //numero bloques capturabilidad MPH
 init_vector yqbloquesmph(1,nqbloquesmph)//yr inicio capturabilidad MPH

//--------------------------------------
// FASES DE ESTIMACI�N
//--------------------------------------
// Capturabilidad
 init_int    opt_qf  //capturabilidad flota
 init_int    opt_qfa  //capturabilidad flota
 init_int    opt_qc  //capturabilidad reclan
 init_int    opt_qmph  //capturabilidad mph
// Selectividad
//(nota:si la selectividad del crucero no se estima se asume igual a 1)
 init_int    optSf_fase  //selectividad flota

 init_int    optS1_fase  //selectividad reclan

// Mortalidad por pesca
 init_int    opt_F

// Desvios de los reclutamientos y No
 init_int    opt_devRt
 init_int    opt_devNo //(dev_No<0 equilibrio)

 init_int    opt_Reclu    //Opci�n BevertonyHolt (1) y Ricker(2)
 init_int    opt_Ro       //Opcion para estimar o dejar fijo log_Ro
 init_number log_priorRo  //Valor fijo para log_Ro, perfiles de verosimilitud

// Puntos biol�gicos de referencia
 init_int    opt_Fpbr //fase de estimaci�n Fpbr (se sugiere sea la �ltima)
 init_int    npbr     //n�mero de pbrs a calcular
 init_vector pbr(1,npbr) //tasas de BDPR ejemplo 30%,40%,60% etc

// Simulacion_proyeccion
 init_int ntime_sim //Años a simular en el futuro
 init_int pR
 init_int oprec
 init_number opt_sim //opcion de simular o estimar (0=simula, 1=estima)
 int reporte_mcmc

 init_vector vec_var(1,nedades)
 init_vector tmp_f(1,nedades)

 init_int test

//#########################################################
INITIALIZATION_SECTION
//#########################################################
  log_Ro           9
  log_Linf         log_Linfprior
  log_k            log_kprior
  log_Lo           log_Loprior
  log_alfa         log_alfa_prior
  log_beta         log_beta_prior
  log_L50f         log_L50priorf
  log_rangof       log_spriorf

  log_L50c         log_L50priorc
  log_rangoc       log_spriorc
  log_M            log_M_prior
  log_qflo         log_prior_qf
  log_qfloa        log_prior_qfa

  log_qcru         log_prior_qc
  log_qmph         log_prior_qmph

//#############################################################
PARAMETER_SECTION
//#############################################################

// Selectividades
//-------------------------------------------------------------
 init_bounded_vector log_L50f(1,nbloques1,-2,2,optSf_fase)
 init_vector log_rangof(1,nbloques1,optSf_fase)

 init_vector log_L50c(1,nbloques2,optS1_fase)
 init_vector log_rangoc(1,nbloques2,optS1_fase)

//-------------------------------------------------------------
// parametros reclutamientos y mortalidades)
//-------------------------------------------------------------
 init_number log_Ro(opt_Ro)
 init_bounded_dev_vector dev_log_Ro(1,ntime,-10,5,opt_devRt)
 init_bounded_vector dev_log_No(1,nedades,-10,10,opt_devNo)
 init_bounded_vector log_F(1,ntime,-20,1.36,opt_F)

//--------------------------------------------------------------
// capturabilidades
//--------------------------------------------------------------
 init_vector log_qflo(1,nqbloques,opt_qf)//cpue industrial
 init_vector log_qfloa(1,nqbloques,opt_qfa)//cpue artesanal
 init_vector log_qcru(1,nqbloquesc,opt_qc)//reclan
 init_vector log_qmph(1,nqbloquesmph,opt_qmph)//mph
//--------------------------------------------------------------
// Crecimiento y mortalidad natural
//--------------------------------------------------------------
 init_number log_Linf(opt_Linf)
 init_number log_k(opt_k)
 init_number log_Lo(opt_Lo)
 init_number log_M(opt_M)
 init_number log_alfa(opt_alfa)
 init_number log_beta(opt_beta)
//--------------------------------------------------------------
// Fpbr
//--------------------------------------------------------------
 init_vector log_Fref(1,npbr,opt_Fpbr);

//*********************************************************************************
// VARIABLES DE ESTADO
//*********************************************************************************
//PRELIMINARY_CALCS_SECTION
//-----------------------------
 vector yrs(1,ntime)
 vector Desemb(1,ntime);
 vector CPUE(1,ntime);
 vector CPUEA(1,ntime);
 vector Bcru(1,ntime);
 vector mph(1,ntime);
 matrix cv_index(1,5,1,ntime);
 matrix nm(1,2,1,ntime);
 vector Unos_edad(1,nedades);
 vector Unos_anos(1,ntime);
 vector Unos_tallas(1,ntallas);
//------------------------------
//FUNCTION Eval_prob_talla_edad
//------------------------------
 number Linf
 number k
 number Lo
 vector mu_edad(1,nedades)
 vector sigma_edad(1,nedades)
 matrix Prob_talla(1,nedades,1,ntallas)
//--------------------------------------
//FUNCTION Eval_selectividad
//--------------------------------------
 matrix S1(1,nbloques1,1,nedades)

 matrix Sel_f(1,ntime,1,nedades)

 matrix S2(1,nbloques2,1,nedades)
 matrix Sel_c(1,ntime,1,nedades);
//--------------------------------------
//FUNCTION Eval_mortalidades
//--------------------------------------
 number M
 matrix F(1,ntime,1,nedades)
 matrix Z(1,ntime,1,nedades)
 matrix S(1,ntime,1,nedades)
//--------------------------------------
//FUNCTION Eval_abundancia
//--------------------------------------
 vector No(1,nedades)
 sdreport_number SSBo
 number phi
 number h
 number alfa
 number beta
 vector Neq(1,nedades);
 matrix N(1,ntime,1,nedades)
 sdreport_vector BD(1,ntime)
 vector Rpred(1,ntime);
 sdreport_vector Reclutas(1,ntime);
 sdreport_vector PHI(1,ntime-1);
//--------------------------------------
//FUNCTION Eval_deinteres
//--------------------------------------
 vector Lf_obs(1,ntime)
 vector Lf_pred(1,ntime)

 vector Lc_obs(1,ntime)
 vector Lc_pred(1,ntime)
 matrix Nv(1,ntime,1,nedades)
 matrix NDv(1,ntime,1,ntallas)
 vector BDo(1,ntime);
 sdreport_vector RPR(1,ntime) //
 sdreport_vector RPRlp(1,ntime)//
//-------------------------------------
//FUNCTION Eval_biomasas
//-------------------------------------
 matrix NMD(1,ntime,1,ntallas)
 matrix NVflo(1,ntime,1,ntallas)
 matrix NVcru(1,ntime,1,ntallas)
 vector BMflo(1,ntime)
 vector BMcru(1,ntime)
 sdreport_vector BT(1,ntime)
 sdreport_vector BV(1,ntime)
//-------------------------------------
//FUNCTION Eval_capturas_predichas
//-------------------------------------
 matrix pred_Ctot_a(1,ntime,1,nedades)
 matrix pred_Ctot(1,ntime,1,ntallas)
 vector pred_Desemb(1,ntime);
 matrix pobs_f(1,ntime,1,ntallas)
 matrix ppred_f(1,ntime,1,ntallas)
 matrix pobsc(1,ntime,1,ntallas)
 matrix ppredc(1,ntime,1,ntallas)
//-------------------------------------
//FUNCTION Eval_indices
//-------------------------------------
 vector pred_CPUE(1,ntime);
 vector pred_CPUEA(1,ntime);
 vector pred_Bcru(1,ntime);
 vector pred_mph(1,ntime);
//-------------------------------------
//FUNCTION Eval_PBR
//-------------------------------------
 vector Fspr(1,nedades)
 vector Zspr(1,nedades)
 vector Nspro(1,nedades)
 vector Nspr(1,nedades)
 vector Nmed(1,nedades)
 number Bspro
 number Bspr
 vector ratio_spr(1,npbr)
 vector Zmed(1,nedades)
 number Bsprmed
 number ratio_Fmed
 number Bo
 vector Brms(1,npbr)
 sdreport_vector RPRrms(1,ntime)
 sdreport_vector Frpr(1,ntime)
//-------------------------------------
//FUNCTION Eval_logverosim
//-------------------------------------
 number suma1
 number suma2
 number suma3
 number penalty
//-------------------------------------
//FUNCTION Eval_funcion_objetivo
//-------------------------------------
 number suma4
 number suma5
 number suma6
 vector likeval(1,10)
//-------------------------------------
//FUNCTION Eval_CTP
//-------------------------------------
 vector Fpbr(1,nedades)
 vector Zpbr(1,nedades)
 vector Np(1,nedades)
 vector Sp(1,nedades)
 number Bp
 number BDp
 vector NMDp(1,ntallas)
 vector CTP(1,ntallas)
 vector CTPp(1,ntallas)
 vector YTP(1,npbr)
 matrix YTPp(1,ntime_sim,1,npbr)
 sdreport_matrix SSBp(1,ntime_sim,1,npbr)// Biomasa desovante proyectada
 matrix BTp(1,ntime_sim,1,npbr)
 number Npplus
 vector Nvp(1,nedades)
 number Nvplus
 vector SDvp(1,ntime_sim)
 vector RPRp(1,npbr)
 matrix Depl(1,ntime_sim,1,npbr)
 vector RPR_last(1,npbr)
 sdreport_vector CBA(1,npbr)
 sdreport_vector CBAp(1,npbr)

 number nm1
 number cuenta1
 number nm2
 number cuenta2

//-------------------------------------

 objective_function_value f

//#########################################################

PRELIMINARY_CALCS_SECTION

//#########################################################

 yrs     = column(data,1);
 Desemb  = column(data,2);
 CPUE    = column(data,4);//cpue industrial
 CPUEA   = column(data,6);
 Bcru    = column(data,8);//reclan
 mph     = column(data,10);//mph
 cv_index(1) = column(data,3);
 cv_index(2) = column(data,5);
 cv_index(3) = column(data,7);
 cv_index(4) = column(data,9);
 cv_index(5) = column(data,11);
 nm(1)       = column(data,12);
 nm(2)       = column(data,13);
 Unos_edad   = 1;// lo uso en operaciones matriciales con la edad
 Unos_anos   = 1;// lo uso en operaciones matriciales con el a�o
 Unos_tallas = 1;// lo uso en operaciones matriciales con el a�o

//#########################################################

RUNTIME_SECTION
  maximum_function_evaluations 5000, 10000, 100000, 500000
  convergence_criteria 1e-6,1e-7,1e-8, 1e-8


//########################################################################

PROCEDURE_SECTION

//########################################################################
 Eval_prob_talla_edad();
 Eval_selectividad();
 Eval_mortalidades();
 Eval_abundancia();
 Eval_biomasas();
 Eval_capturas_predichas();
 Eval_indices();
 Eval_PBR();
 Eval_deinteres();
 Eval_logverosim();
 Eval_funcion_objetivo();

 if(last_phase()){
 	Eval_CTP();
 }

  if(mceval_phase())
    {
    ofstream mcmc_report("mcmc2.mcmc",ios::app);
    mcmc_report << BD(ntime) << " " << BV(ntime) << endl;
    mcmc_report.close();

    ofstream doris("ctp.mcmc",ios::app);
    doris << CBAp(1) << " " << CBAp(2) << " " << CBAp(3) << endl;
    doris.close();

    }


//========================================================================

FUNCTION Eval_prob_talla_edad

//========================================================================
 int i, j;
// genero una clave edad-talla para otros calculos. Se modela desde L(1)
//---------------------------------------------------------
 Linf = mfexp(log_Linfprior);
 k    = mfexp(log_kprior);
 Lo   = mfexp(log_Lo);

   mu_edad(1) = Lo;
   for (i=2;i<=nedades;i++){
        mu_edad(i) = mu_edad(i-1)*mfexp(-k)+Linf*(1-mfexp(-k)); //Shnute y Fornier 1980
         }
   sigma_edad = mfexp(log_alfa)+mfexp(log_beta)*mu_edad;
   //sigma_edad = mfexp(log_beta)*mu_edad;
   //sigma_edad = vec_var;
   Prob_talla = ALK( mu_edad, sigma_edad, Tallas);


//========================================================================

FUNCTION dvar_matrix ALK(dvar_vector& mu, dvar_vector& sig, dvector& x)

//========================================================================
	//RETURN_ARRAYS_INCREMENT();
	int i, j;
	dvariable z1;
	dvariable z2;
	int si,ni; si=mu.indexmin(); ni=mu.indexmax();
	int sj,nj; sj=x.indexmin(); nj=x.indexmax();
	dvar_matrix pdf(si,ni,sj,nj);
	pdf.initialize();
	double xs=0.5*(x[sj+1]-x[sj]);
	for(i=si;i<=ni;i++) //loop over ages
	{
		 for(j=sj;j<=nj;j++) //loop over length bins
		{
			z1=((x(j)-xs)-mu(i))/sig(i);
			z2=((x(j)+xs)-mu(i))/sig(i);
			pdf(i,j)=cumd_norm(z2)-cumd_norm(z1);
		}//end nbins
		pdf(i)/=sum(pdf(i));
	}//end nage
	//RETURN_ARRAYS_DECREMENT();
	return(pdf);

//========================================================================

FUNCTION Eval_selectividad

//========================================================================
 int i,j;
//***********************************************************************
// SELECTIVIDAD LOGISTICA
//***********************************************************************
//-----------------------------------------------------------------------
// FLOTA
//-----------------------------------------------------------------------
 for (j=1;j<=nbloques1;j++){
   S1(j)=1/(1+exp(-log(19)*(edades-exp(log_L50f(j)))/exp(log_rangof(j))));}

   for (i=1;i<=ntime;i++){
      for (j=1;j<=nbloques1;j++){
            if (yrs(i)>=ybloques1(j)){
                  Sel_f(i)=S1(j);}}}

	//for (i=1;i<=ntime;i++){
	//    Sel_f(i) = tmp_f;}


//-----------------------------------------------------------------------
// CRUCEROS RECLAN
//-----------------------------------------------------------------------
    Sel_c=1.0;

  if(active(log_L50c)){

  for (j=1;j<=nbloques2;j++){
     S2(j)=1/(1+exp(-log(19)*(edades-exp(log_L50c(j)))/exp(log_rangoc(j))));}

   for (i=1;i<=ntime;i++){
      for (j=1;j<=nbloques2;j++){
           if (yrs(i)>=ybloques2(j)){
                Sel_c(i)=S2(j);}}}
  }


//========================================================================

FUNCTION Eval_mortalidades

//========================================================================

 M=exp(log_M);
 F=elem_prod(Sel_f,outer_prod(mfexp(log_F),Unos_edad));
 Z=F+M;
 S=mfexp(-1.0*Z);




//========================================================================
FUNCTION Eval_abundancia
//========================================================================
 int i, j;

  if(opt_Ro<0)
  {
  log_Ro=log_priorRo;
  }

 // Biomasa desovante virgen de largo plazo
  No(1)=mfexp(log_Ro);  //  +0.5*square(cvar(1)));
 for (int j=2;j<=nedades;j++)
 {No(j)=No(j-1)*mfexp(-1.*M);}
  No(nedades)=No(nedades)/(1-exp(-1.*M));
  //No(nedades)=No(nedades)*exp(-1.*M)/(1-exp(-1.*M)); //grupo plus
  SSBo=sum(elem_prod(No*exp(-dt(1)*M)*Prob_talla,elem_prod(msex,Wmed)));
  phi=SSBo/mfexp(log_Ro);


//Stock-Recluta
 h=mfexp(log_h_prior);
 if(opt_Reclu==1){
 alfa=4*h*mfexp(log_Ro)/(5*h-1);//
 beta=(1-h)*SSBo/(5*h-1);
 }else{
 alfa=1.25*log(5*h)-log(phi);
 beta=1.25*log(5*h)/SSBo;}

// genero una estructura inicial en equilibrio para el primer a�o
 //Neq(1)=mfexp(log_Ro+0.5*square(cvar(1)));
    Neq(1)=No(1);
 for (j=2;j<=nedades;j++){
   Neq(j)=Neq(j-1)*mfexp(-Z(1,j-1));}
   Neq(nedades)=Neq(nedades)/(1-mfexp(-1.*Z(1,nedades)));

 // Abundancia inicial
 N(1)     = elem_prod(Neq,mfexp(dev_log_No));
 BD(1)    = sum(elem_prod(elem_prod(N(1),mfexp(-dt(1)*Z(1)))*Prob_talla,elem_prod(msex,Wmed)));
 Rpred(1) = N(1,1);//

// se estima la sobrevivencia por edad(a+1) y a�o(t+1)
 for (i=1;i<ntime;i++)
 {
     if(opt_Reclu==1){
     Rpred(i+1)=alfa*BD(i)/(beta+BD(i));
	 }
      else if(opt_Reclu==2) {
		 Rpred(i+1)=BD(i)*mfexp(alfa-1.*beta*BD(i));
	 } else {
		 Rpred(i+1)= exp(log_Ro);
	 }

     N(i+1,1)=Rpred(i+1)*mfexp(dev_log_Ro(i+1));// Reclutas
     N(i+1)(2,nedades)=++elem_prod(N(i)(1,nedades-1),S(i)(1,nedades-1));
     N(i+1,nedades) = N(i+1,nedades)+N(i,nedades)*S(i,nedades);// grupo plus
     BD(i+1)=sum(elem_prod(elem_prod(N(i+1),exp(-dt(1)*Z(i+1)))*Prob_talla,elem_prod(msex,Wmed)));

	 PHI(i) = BD(i)/N(i+1,1);
 }
     Reclutas=column(N,1);
     //PHI=elem_div(BD,Reclutas);

//========================================================================

FUNCTION Eval_biomasas

//========================================================================

// NMD = elem_prod(N,mfexp(-dt(1)*Z))*Prob_talla;
 //NMD = elem_prod(NMD,outer_prod(Unos_anos,msex));

 NVflo   = elem_prod(elem_prod(N,mfexp(-Z)),Sel_f)*Prob_talla;          //FLOTA
 NVcru   = elem_prod(elem_prod(N,mfexp(-dt(2)*(Z))),Sel_c)*Prob_talla;  //RECLAN

// vectores de biomasas derivadas
 BMflo   = NVflo*Wmed;
 BMcru   = NVcru*Wmed;
 BV      = BMflo;
 BT      = N*Prob_talla*Wmed;


//========================================================================

FUNCTION Eval_capturas_predichas

//========================================================================

// matrices de capturas predichas por edad y a�o
 pred_Ctot_a = elem_prod(elem_div(F,Z),elem_prod(1.-S,N));
 pred_Ctot   = pred_Ctot_a*Prob_talla;
 pred_Desemb = pred_Ctot*Wmed;

// PROPORCIONES  matrices de proporcion de capturas por talla y a�o
 pobs_f      = elem_div(Ctot,outer_prod(rowsum(Ctot+1e-10),Unos_tallas));
 ppred_f     = elem_div(pred_Ctot,outer_prod(rowsum(pred_Ctot+1e-10),Unos_tallas));
 pobsc       = elem_div(Ncru,outer_prod(rowsum(Ncru+1e-10),Unos_tallas));
 ppredc      = elem_div(NVcru,outer_prod(rowsum(NVcru+1e-10),Unos_tallas));

//========================================================================

FUNCTION Eval_indices

//========================================================================
    for (int i=1;i<=ntime;i++){
      for (int j=1;j<=nqbloques;j++){
              if (yrs(i)>=yqbloques(j)){
                 pred_CPUE(i)=mfexp(log_qflo(j))*BMflo(i);}
       }
   }

    for (int i=1; i<=ntime; i++){
	  for (int j=1; j<=nqabloques; j++){
		      if (yrs(i)>=yqabloques(j)){
				  pred_CPUEA(i)=mfexp(log_qfloa(j))*BMflo(i);
										}
									   }
								}

   for (int i=1;i<=ntime;i++){
      for (int j=1;j<=nqbloquesc;j++){
              if (yrs(i)>=yqbloquesc(j)){
                 pred_Bcru(i)=mfexp(log_qcru(j))*BMcru(i);}
       }
   }


   for (int i=1;i<=ntime;i++){
      for (int j=1;j<=nqbloquesmph;j++){
              if (yrs(i)>=yqbloquesmph(j)){
                 pred_mph(i)=mfexp(log_qmph(j))*BD(i);}
       }
   }

//========================================================================

FUNCTION Eval_PBR

//========================================================================

 if(opt_Ro<0)
 {
        log_Ro=log_priorRo;
 }

 for (int i=1;i<=npbr;i++){

 Fspr=Sel_f(ntime)*mfexp(log_Fref(i));
 Zspr=Fspr+M;

 Nspro(1)=mfexp(log_Ro);  //+0.5*square(cvar(1)));
 Nspr(1)=mfexp(log_Ro);  //+0.5*square(cvar(1)));

 for (int j=2;j<=nedades;j++)
 {
   Nspro(j)=Nspro(j-1)*mfexp(-1.*M);
   Nspr(j)=Nspr(j-1)*mfexp(-Zspr(j-1));
 }

   Nspro(nedades)=Nspro(nedades)/(1-mfexp(-1.*M));
   Nspr(nedades)=Nspr(nedades)/(1-mfexp(-Zspr(nedades)));

 Bspro   = sum(elem_prod(Nspro*mfexp(-dt(1)*M)*Prob_talla,elem_prod(msex,Wmed)));
 Bspr    = sum(elem_prod(elem_prod(Nspr,mfexp(-dt(1)*Zspr))*Prob_talla,elem_prod(msex,Wmed)));

 //BD_lp=sum(elem_prod(elem_prod(Neq,mfexp(-dt(1)*Zpbr))*Prob_talla,elem_prod(msex,Wmed)));

 ratio_spr(i)=Bspr/Bspro;
 Bo = Bspro;
 Brms(i)=Bo*(ratio_spr(i)-0.05);
 }
 RPRrms = BD/Brms(1);
 Frpr   = exp(log_F)/mfexp((log_Fref(1)));

//========================================================================

FUNCTION Eval_deinteres

//========================================================================

//------------------------------------------
// TALLAS PROMEDIO
//------------------------------------------
 Lf_obs   = Tallas*trans(pobs_f);
 Lf_pred  = Tallas*trans(ppred_f);
 Lc_obs   = Tallas*trans(pobsc);
 Lc_pred  = Tallas*trans(ppredc);


//------------------------------------------
// Rutina para calcular RPR
//------------------------------------------

 Nv=N;// solo para empezar los calculos

// se estima la sobrevivencia por edad(a+1) y ano(t+1)
 for (int i=1;i<ntime;i++)
 {
     Nv(i+1)(2,nedades)=++Nv(i)(1,nedades-1)*exp(-1.0*M);
    // Nv(i+1,nedades)=Nv(i+1,nedades)+Nv(i,nedades)*exp(-1.0*M);// grupo plus
 }

 NDv   = elem_prod((Nv*exp(-dt(1)*M))*Prob_talla,outer_prod(Unos_anos,msex));
 BDo   = NDv*Wmed;
 RPR   = elem_div(BD,BDo);
 RPRlp = BD/SSBo;

 //========================================================================

FUNCTION Eval_logverosim

//========================================================================

// esta funcion evalua el nucleo de las -log-verosimilitudes marginales para
// series con datos 0.
 int i;

 suma1=0; suma2=0; suma3=0; suma4=0; penalty=0;

 for (i=1;i<=ntime;i++)
 {
   if (CPUE(i)>0){
    suma1+=square(log(CPUE(i)/pred_CPUE(i))*1/cv_index(2,i));}
	if (CPUEA(i)>0){
    suma2+=square(log(CPUEA(i)/pred_CPUEA(i))*1/cv_index(3,i));}
  if (Bcru(i)>0){
    suma3+=square(log(Bcru(i)/pred_Bcru(i))*1/cv_index(4,i));}
  if (mph(i)>0){
    suma4+=square(log(mph(i)/pred_mph(i))*1/0.1);}
  }

//========================================================================

FUNCTION Eval_funcion_objetivo

//========================================================================

 suma5=0; suma6=0; penalty=0;
 int i;

 likeval(1)=0.5*suma1;//CPUE Industrial
 likeval(2)=0.5*suma2;//CPUE Artesanal
 likeval(3)=0.5*suma3;//RECLAN
 likeval(4)=0.5*suma4;//MPH
 likeval(5) = 0.5*norm2(elem_div(log(elem_div(Desemb,pred_Desemb)),cv_index(1))); //Desemb

 suma5 = sum(nm(1)*elem_prod(pobs_f,log(ppred_f)));
 suma6 = sum(nm(2)*elem_prod(pobsc,log(ppredc)));

 likeval(6) = -1.*suma5; //FLOTA
 likeval(7) = -1.*suma6; //RECLAN


// lognormal Ninicial y Reclutas
 if(active(dev_log_Ro)){
 likeval(8)=1./(2*square(cvar(1)))*norm2(dev_log_Ro);}

 if(active(dev_log_No)){
 likeval(9)=1./(2*square(cvar(2)))*norm2(dev_log_No);}

 if (active(log_Lo)){
 likeval(10)=0.5*square((log_Loprior-log_Lo)/cv_priors(3));}


//--------------------------------------------------------------
//PENALIZACIONES
//--------------------------------------------------------------
 if(active(log_M)){
 penalty+=100*(square(log_M_prior-log_M));}

 //if(active(log_qflo)){
 //penalty+=0.5*norm2((log_qflo-log_prior_qf)/cv_qf);}
 //if(active(log_qfloa)){
 //penalty+=0.5*norm2((log_qfloa-log_prior_qfa)/cv_qfa);}
 //if(active(log_qcru)){
 //penalty+=0.5*norm2((log_qcru-log_prior_qc)/cv_qcru);}
 // if(active(log_qmph)){
 // penalty+=0.5*norm2((log_qmph-log_prior_qmph)/cv_qmph);}


 if(active(log_Fref)){
 penalty+=1000*norm2(ratio_spr-pbr);}

 f=opt_sim*(sum(likeval)+penalty);

  //depuro(penalty)

//========================================================================

FUNCTION Eval_CTP

//========================================================================

// re-calcula la CBA para el �ltimo a�o dado los PBR entregados

 for (int j=1;j<=npbr;j++){
     Fpbr   = Sel_f(ntime)*exp(log_Fref(j));//
     Zpbr   = Fpbr+M;
     CTP    = elem_prod(elem_div(Fpbr,Zpbr),elem_prod(1.-exp(-1.*Zpbr),N(ntime)))*Prob_talla;
     YTP(j) = sum(elem_prod(CTP,Wmed));
     }

 for (int j=1;j<=npbr;j++){

 Np   = N(ntime);
 Sp   = S(ntime);
 BDp  = BD(ntime);
 Fpbr = F(ntime);//
 Zpbr = Z(ntime);

   for (int i=1;i<=ntime_sim;i++){

//*****************************************************************************
// ESTIMACION DEL RECLUTAMIENTO Ejemplos
//*****************************************************************************
   if(oprec==1){//reclutamiento para proyeccion segun relacion stock-recluta o no
     if(opt_Reclu==1){ //Beverton Y Holt
       Np(1)=pR*alfa*BDp/(beta+BDp)*mfexp(0.5*square(cvar(1)));}
     if(opt_Reclu==2){ //Ricker
       Np(1)=pR*BDp*mfexp(alfa-1.0*beta*BDp)*mfexp(0.5*square(cvar(1)));}
     if(opt_Reclu==3){ //sin relacion stock-recluta
       Np(1)=pR*mfexp(log_Ro);}}  //+0.5*square(cvar(1)))

 //Periodo de baja productividad adecuados para anchoveta centro-norte
   if(oprec==2){//reclutamiento proyectado promedio periodo x
      Np(1)=mean(Reclutas(ntime-5,ntime-1));}
 //Periodo de alta productividad
   if(oprec==3){//reclutamiento proyectado promedio hasta 1998
      Np(1)=mean(Reclutas(1,ntime-18));}
//*****************************************************************************

  Np(2,nedades)=++elem_prod(Np(1,nedades-1),Sp(1,nedades-1));
  Np(nedades)+=Np(nedades)*Sp(nedades);

  Fpbr = Sel_f(ntime)*exp(log_Fref(j));//
  Zpbr = Fpbr+M;
  Sp=exp(-1.*Zpbr);

  Bp=sum(elem_prod(Np*Prob_talla,Wmed));

  NMDp = elem_prod(Np,mfexp(-dt(1)*(Zpbr)))*Prob_talla;
  BDp  = sum(elem_prod(elem_prod(NMDp,msex),Wmed));
  CTPp = elem_prod(elem_div(Fpbr,Zpbr),elem_prod(Np,(1-Sp)))*Prob_talla;

  YTPp(i,j)=sum(elem_prod(CTPp,Wmed));
  SSBp(i,j)=BDp;
  BTp(i,j)=Bp;

 }

// proyectado desde RECLAN para el mismo ano
  CBA(j) = YTP(j);

// proyectado desde fin de ano al proximo
  CBAp(j) = YTPp(1,j);// es para el ano proyectado revisar si es 1 o 2

 }

 for (int i=1;i<=ntime_sim;i++)
  {
     Nvplus=Nvp(nedades)*exp(-1.0*M);
     Nvp(2,nedades)=++Nvp(1,nedades-1)*exp(-1.0*M);
     Nvp(nedades)+=Nvplus;
     Nvp(1)= exp(log_Ro);
     SDvp(i)=sum(elem_prod(Nvp*Prob_talla,elem_prod(Wmed,msex)));
  }


//#########################################################

REPORT_SECTION

//#########################################################

 report << "YRS" << endl;
 report << yrs << endl;
//----------------------------------------
// INDICES OBSERVADOS VS PREDICHOS
//----------------------------------------
 report << "desemb" << endl;
 report << Desemb << endl;
 report << "desemb_pred" << endl;
 report << pred_Desemb << endl;
 report << "reclan" << endl;
 report << Bcru << endl;
 report << "reclan_pred" << endl;
 report << pred_Bcru << endl;
 report << "CPUE" << endl;
 report << CPUE << endl;
 report << "pred_cpue" << endl;
 report << pred_CPUE << endl;
 report << "CPUEA" << endl;
 report << CPUEA << endl;
 report << "pred_cpuea" << endl;
 report << pred_CPUEA << endl;
 report << "mph" << endl;
 report << mph << endl;
 report << "mph_pred" << endl;
 report << pred_mph << endl;
//-----------------------------------------
// INDICADORES POBLACIONALES
//-----------------------------------------
 report << "BD" << endl;
 report << BD << endl;
 report << "BT" << endl;
 report << BT << endl;
 report << "NV" << endl;
 report << NVflo << endl;
 report << "BV" << endl;
 report << BMflo << endl;
 report << "Reclutas" << endl;
 report << Rpred << endl;
 report << column(N,1)<< endl;
 report << "PHI" << endl;
 report << PHI << endl;
 report << "pred_Ctot" << endl;
 report << pred_Ctot << endl;

 report << "F" << endl;
 report << rowsum(F)/nedades<< endl;
//------------------------------------------
// TALLAS PROMEDIO
//------------------------------------------
 report << "Lf_obs" << endl;
 report << Lf_obs<< endl;
 report << "Lf_pred" << endl;
 report << Lf_pred<< endl;
 report << "Lc_obs" << endl;
 report << Lc_obs<< endl;
 report << "Lc_pred" << endl;
 report << Lc_pred<< endl;
 //-------------------------------------------
// SELECTIVIDADES
//-------------------------------------------
 report << "Sflo_age" <<endl;
 report << Sel_f << endl;
 report << "Scru_age" << endl;
 report << Sel_c << endl;
 //-------------------------------------------
// PROPORCION DE LAS CAPTURAS
//-------------------------------------------
 report << "pf_obs" << endl;
 report << (pobs_f)<< endl;
 report << "pf_pred" << endl;
 report << (ppred_f)<< endl;
 report << "pcru_obs" << endl;
 report << (pobsc)<< endl;
 report << "pcru_pred" << endl;
 report << (ppredc)<< endl;
 //------------------------------------------
 report << "N" << endl;
 report << N << endl;
 report << "Capt_age" <<endl;
 report << pred_Ctot_a <<endl;
 report << "F_age" <<endl;
 report << F << endl;
//-----------------------------------------
// INDICADORES DE REDUCCI�N DEL STOCK
//-----------------------------------------
 report << "BDo" << endl;
 report << BDo << endl;
 report << "SSBo" << endl;
 report << SSBo << endl;
 report << "RPR" << endl;
 report << RPR << endl;
 report << "RPRlp" << endl;
 report << RPRlp << endl;
 report << "RPRrms"<<endl;
 report <<  RPRrms << endl;
 report << "Frpr"<<endl;
 report <<  Frpr <<endl;
 report << "log_Ro" <<endl;
 report << log_Ro <<endl;
//------------------------------------------
// CRECIMIENTO Y MORTALIDAD NATURAL
//------------------------------------------
 report << "Lo"  << endl;
 report << exp(log_Lo) << endl;
 report << "M" << endl;
 report << M <<endl;
 report << "mu_edad" << endl;
 report << mu_edad << endl;
 report << "Prob_talla" << endl;
 report << Prob_talla << endl;
//------------------------------------------
// VEROSIMILITUD
//------------------------------------------
 report << "Like" <<endl;
 report << likeval << endl;
//------------------------------------------
// CAPTURABILIDADES
//------------------------------------------
 report << "q_flo" <<endl;
 report << exp(log_qflo) << endl;
 report << "q_floa" <<endl;
 report << exp(log_qfloa) << endl;
 report << "q_cru_reclan" <<endl;
 report << exp(log_qcru) << endl;

//------------------------------------------
//PARAMETROS RELACI�N STOCK-RECLUTA
//-----------------------------------------
 report << "alfa" << endl;
 report << alfa << endl;
 report << "beta" <<endl;
 report << beta << endl;
 report << "phi" << endl;
 report << phi  << endl;
//----------------------------------------
// PUNTOS BIOL�GICOS DE REFERENCIA
//----------------------------------------
 report << "pSPR_Fpbr" <<endl;
 report << ratio_spr << endl;
 report << "Frms"<<endl;
 report << mfexp((log_Fref(1))) << endl;
 report << "Bo"<<endl;
 report << Bo << endl;
 report << "Brms" << endl;
 report << Brms << endl;
 report << "msex_edad"<<endl;
 report << msex << endl;
 report << "Wmed_edad"<<endl;
 report << Wmed << endl;
 report << "Nspr"<<endl;
 report << Nspr << endl;
 report << "Nspro"<<endl;
 report << Nspro << endl;
//-------------------------------------
// CBA Y PROYECCIONES
//--------------------------------------
 report << "RPR_lp" << endl;
 report << RPRp << endl;
 report << "BT_proy" << endl;
 report << BTp << endl;
 report << "BD_proy" << endl;
 report << SSBp << endl;
 report << "CBA" << endl;
 report << CBA<< endl;
 report << "CBAp" << endl;
 report << CBAp << endl;
 report << "Np" << endl;
 report << Np  << endl;
 report << "Nvp" << endl;
 report << CTPp  << endl;
 report << "CTP" << endl;
 report << CTP  << endl;


//########################################################################

FINAL_SECTION

//########################################################################

 time(&finish);
 elapsed_time=difftime(finish,start);
 hour=long(elapsed_time)/3600;
 minute=long(elapsed_time)%3600/60;
 second=(long(elapsed_time)%3600)%60;
 cout<<endl<<endl<<"*********************************************"<<endl;
 cout<<"--Start time:  "<<ctime(&start)<<endl;
 cout<<"--Finish time: "<<ctime(&finish)<<endl;
 cout<<"--Runtime: ";
 cout<<hour<<" hours, "<<minute<<" minutes, "<<second<<" seconds"<<endl;
 cout<<"*********************************************"<<endl;

 //######################################################################
