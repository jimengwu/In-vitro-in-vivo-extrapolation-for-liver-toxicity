#---------------------------13nm code----------
humanPBPK.code <- '
$PARAM @annotated
// -------------------1. Physiological parameters for weight of organs---------
Body_weight : 73 : kg, Bodyweight https://www.sciencedirect.com/science/article/pii/S0169409X01001508?via%3Dihub

F_W_Liver  : 0.0257  : (Unitless, Fractional liver tissue (Brown 1997) table 7, https://journals.sagepub.com/doi/epdf/10.1177/074823379701300401)
F_W_Brain  : 0.02  : Unitless, Fractional brain tissue (Brown 1997) table 7
F_W_Lung   : 0.0076  : Unitless, Fractional lung tissue (Brown 1997) table 7
F_W_Kidney : 0.0044  : Unitless, Fractional kidney tissue (Brown 1997) table 7
F_W_Spleen : 0.0026  : Unitless, Fractional spleen tissue (Davies 1993) table 7
F_W_GI : 0.0171  : (Unitless, Fractional GI tract tissue (Brown 1997) table 7
F_W_Plasma : 0.073  : https://doi.org/10.1016/S0146-6453(03)00002-2 5.3kg of blood, 73kg of body weight
//F_W_Plasma : 0.029  : Unitless, Fractional plasma (Davies 1993)
// -------------------1. Physiological parameters for weight of organs---------

// -------------------2. Physiological parameters for blood flow of organs---------
// QC : 13.98  : mL/min, Cardiac output 
QC : 5200 * 60 / 1000  : L/h https://doi.org/10.1177/074823379701300401  table22
F_Q_Brain : 0.114  : Unitless, Fraction blood flow to brain of the cardiac output (Brown 1997) table 23
F_Q_Liver : 0.227  : Unitless, Fraction blood flow to liver of the cardiac output (Brown 1997) table 23
F_Q_Kidney : 0.175  : Unitless, Fraction blood flow to kidney (Brown 1997) table 23
F_Q_Spleen : 0.0138  : https://link.springer.com/article/10.1023/a:1018943613122
F_Q_GI : 0.1387  :https://link.springer.com/article/10.1007/s10928-011-9232-2#ref-CR13
// -------------------2. Physiological parameters for blood flow of organs---------

// ------------------3. Partition coefficient between tissue and blood---------------
P_Liver  : 0.08   :unitless
P_Brain  : 0.15  :unitless
P_Kidney : 0.15  :unitless
P_Spleen : 0.15  :unitless
P_Lung   : 0.15  :unitless
P_Rest   : 0.15  :unitless
P_GI     : 0.15  :unitless
// ------------------3. Partition coefficient between tissue and blood---------------

// -------------------4. Physiological parameters for blood volume fraction in tissues---------------
F_BV_Liver : 0.11  : Unitless, Fraction blood volume in liver, percentage of liver (Brown 1997) table 30
F_BV_GI : 0.029  :https://link.springer.com/article/10.1007/s10928-011-9232-2#ref-CR13
F_BV_Brain : 0.04  : Unitless, Fraction blood volume in brain, percentage of brain (Brown 1997) table 30
F_BV_Kidney : 0.36  : Unitless, Fraction blood volume in kidney, percentage of kidney (Brown 1997) table 30
F_BV_Spleen : 0.0516  : https://link.springer.com/article/10.1007/s10928-011-9232-2#ref-CR13
F_BV_Lung : 0.1  : https://link.springer.com/article/10.1007/s10928-011-9232-2#ref-CR13
F_BV_Rest : 0.076  : https://link.springer.com/article/10.1007/s10928-011-9232-2#ref-CR13
// -------------------4. Physiological parameters for blood volume fraction in tissues---------------


// -------------------5. Diffusion limitation coefficient constants (Permeability coefficient)------------------------
DLC_Liver : 0.001  : unitless
DLC_Brain : 0.000001  : due to blood brain barrier unitless
DLC_Kidney : 0.001 :unitless
DLC_Spleen : 0.03 :unitless
DLC_Lung : 0.001:unitless
DLC_GI : 0.001:unitless
DLC_Rest : 0.000001:unitless
// -------------------5. Diffusion limitation coefficient constants (Permeability coefficient)------------------------


// -------------------6. Endocytosis-related parameters------------------------

K_release_Liver : 0.001  : h-1, 13nm
K_uptake_Liver     :  20    : h-1, 13nm
K_50_Liver      :   48   : h, 13nm
n_Liver         :   5    : Unitless, 13nm


K_release_GI    :    0.001    : h-1
K_uptake_GI        :    0.075    : h-1
K_50_GI         :    24       :h
n_GI            :     5       :Unitless



K_release_Spleen : 0.001  : h-1, 13nm
K_uptake_Spleen     : 40     : h-1, 13nm
K_50_Spleen      : 48     :h, 13nm
n_Spleen         : 5      :Unitless, 13nm



K_release_Kidney : 0.0004 : h-1, 13nm
K_uptake_Kidney     : 0.075  : h-1, 13nm
K_50_Kidney      : 24     :h
n_Kidney         : 5      :Unitless, 13nm

K_release_Lung   : 0.003 : h-1, 13nm
K_uptake_Lung       : 0.075 : h-1, 13nm
K_50_Lung        : 24    :h, 13nm
n_Lung           : 5     :Unitless, 13nm

A_cap_liver : 1000:μg/kg tissue

K_GI_b  :   4e-5     : h-1, absorption rate of GI tract
//A_cap_spleen : 4000:μg/kg tissue
//A_cap_kidney : 4000:μg/kg tissue
// -------------------6. Endocytosis-related parameters------------------------


// -------------------7. Excretion----------------------------
Kbile  : 0.00003   : Biliary clearance (L/hr), 13nm
Kurine : 0.000003  : Urine clearance (L/hr), 13nm
Kfecal : 0.000003  : fitting L/h
// -------------------7. Excretion----------------------------

$MAIN
// -------------------1. Blood flow------------

double Q_Brain = F_Q_Brain * QC  ;// L/h, Blood flow to brain
double Q_Kidney = F_Q_Kidney * QC ; // L/h, Blood flow to kidney
double Q_Liver = F_Q_Liver * QC ; // L/h, Blood flow to liver
double Q_Spleen = F_Q_Spleen * QC ; // L/h, Blood flow to Spleen
double Q_GI = F_Q_GI * QC  ;// L/h, Blood flow to GI

//double Q_Rest =  QC - Q_Kidney - Q_Liver - Q_Brain - Q_Spleen ;// L/h, Blood flow to the rest of body
//double Q_balance =  QC - Q_Kidney - Q_Liver - Q_Brain - Q_Spleen - Q_Rest;// L/h, Blood flow to the rest of body

double Q_Rest =  QC - Q_Kidney - Q_Liver - Q_Brain - Q_Spleen - Q_GI ;// L/h, Blood flow to the rest of body
double Q_balance =  QC - Q_Kidney - Q_Liver - Q_Brain - Q_Spleen - Q_Rest-Q_GI;// L/h, Blood flow to the rest of body


// ------------------1. Blood flow--------------


// --------------2. Tissue Volume-------------
double W_Liver = Body_weight * F_W_Liver;  // Liver, kg
double W_Brain = Body_weight * F_W_Brain;  // Brain
double W_Kidney = Body_weight * F_W_Kidney;  // Kidney
double W_Spleen = Body_weight * F_W_Spleen; // Spleen
double W_Lung = Body_weight * F_W_Lung;  // Lungs
double W_GI = Body_weight * F_W_GI;  // GI

double W_Plasma = Body_weight * F_W_Plasma;  // plasma


//double W_Rest = Body_weight - W_Liver - W_Brain - 
//                W_Kidney - W_Spleen - W_Lung - W_Plasma;

double W_Rest = Body_weight - W_Liver - W_Brain - 
                W_Kidney - W_Spleen - W_Lung - W_Plasma-W_GI;

//double W_balance = Body_weight
//    - W_Liver
//    - W_Brain
//    - W_Kidney
//    - W_Spleen
//    - W_Lung
//    - W_Plasma
//    - W_Rest;

double W_balance = Body_weight
    - W_Liver
    - W_Brain
    - W_Kidney
    - W_Spleen
    - W_Lung
    - W_Plasma
    - W_Rest
    -W_GI;

double W_Liver_b = W_Liver * F_BV_Liver; // Weight/volume of capillary blood in liver compartment
double W_Liver_t = W_Liver - W_Liver_b;  // Weight/volume of tissue in liver compartment

double W_GI_b = W_GI * F_BV_GI;  // Weight/volume of capillary blood in GI compartment
double W_GI_t = W_GI - W_GI_b;  // Weight/volume of tissue in GI compartment


double W_Brain_b = W_Brain * F_BV_Brain; // Weight/volume of capillary blood in Brain compartment
double W_Brain_t = W_Brain - W_Brain_b;  // Weight/volume of tissue in Brain compartment


double W_Kidney_b = W_Kidney * F_BV_Kidney; // Weight/volume of capillary blood in Kidney compartment
double W_Kidney_t = W_Kidney - W_Kidney_b;  // Weight/volume of tissue in Kidney compartment

double W_Spleen_b = W_Spleen * F_BV_Spleen; // Weight/volume of capillary blood in Spleen compartment
double W_Spleen_t = W_Spleen - W_Spleen_b;  // Weight/volume of tissue in Spleen compartment

double W_Lung_b = W_Lung * F_BV_Lung;  // Weight/volume of capillary blood in Lung compartment
double W_Lung_t = W_Lung - W_Lung_b;  // Weight/volume of tissue in Lung compartment

double W_Rest_b = W_Rest * F_BV_Rest;  // Weight/volume of capillary blood in Rest compartment
double W_Rest_t = W_Rest - W_Rest_b;  // Weight/volume of tissue in Rest compartment

// --------------2. Tissue Volume-------------

//----------------3. Endocytosis parameter-------------
//double K_up_Liver = K_uptake_Liver;
double K_up_Spleen = K_uptake_Spleen;
//double K_up_Liver  = (K_uptake_Liver  * pow(TIME,n_Liver))  / (pow(K_50_Liver, n_Liver)    + pow(TIME,n_Liver));
//double K_up_Spleen = (K_uptake_Spleen * pow(TIME,n_Spleen)) / ( pow(K_50_Spleen, n_Spleen) + pow(TIME,n_Spleen));

double K_up_Kidney = K_uptake_Kidney;
double K_up_Lung   = K_uptake_Lung; 
//double K_up_GI   = (K_uptake_GI     * pow(TIME,n_GI))     / ( pow(K_50_GI, n_Lung)       + pow(TIME,n_Lung));
double K_up_GI   = K_uptake_GI;


//----------------3. Endocytosis parameter-------------

$CMT MBA MBV M_Lung_b M_Lung_t M_Lung_pc M_Brain_b M_Brain_t M_Rest_b M_Rest_t 
M_Kidney_b M_Kidney_t M_Kidney_pc M_Spleen_b M_Spleen_t M_Spleen_pc M_Liver_b 
M_Liver_t M_Liver_pc M_GI_b M_GI_pc M_GI_t M_GI_lumen M_Urine M_bile M_fecal 
AUC_Lt AUC_Kt AUC_St AUC_Lut AUC_Plasma

$ODE

//---------------------------DIFFERENTIAL EQUATION-----------------------------
// -----------------------1. Blood: Concentrations in the blood compartment---------------

// --------------------------- 1.1 Arterial Blood part ---------------------------------------
double K_up_Liver  = K_uptake_Liver * (1 - M_Liver_pc/(A_cap_liver  * W_Liver_t));  

//double K_up_Spleen = K_uptake_Spleen * (1 - M_Spleen_pc/(A_cap_spleen  * W_Spleen_t));  
//double K_up_Kidney = K_uptake_Kidney * (1 - M_Kidney_pc/(A_cap_kidney  * W_Kidney_t));  
double RBA = QC * (M_Lung_b / W_Lung_b) - QC * (MBA / (W_Plasma * 0.2));  // rate of change in Arterial blood
dxdt_MBA = RBA; // MBA with unit of mg  
// --------------------------- 1.1 Arterial Blood part ---------------------------------------

// --------------------------- 1.2 Venous Blood part ---------------------------------------


double RBV = (Q_Liver * (M_Liver_b / W_Liver_b) + Q_Brain * (M_Brain_b / W_Brain_b) 
            + Q_Kidney * (M_Kidney_b / W_Kidney_b) + Q_Rest * (M_Rest_b / W_Rest_b)  
            - QC * (MBV / (W_Plasma * 0.8)));
dxdt_MBV = RBV;
double C_Plasma = (MBA+ MBV) / W_Plasma;
dxdt_AUC_Plasma = C_Plasma*1000;
// --------------------------- 1.2 Venous Blood part ---------------------------------------

// --------------------------------- 2. Lung compartment------------------------------------

// ----------------------------------2.1 lung blood compartment ----------------------------
double R_Lung_b = (QC * ((MBV / (W_Plasma * 0.8)) - (M_Lung_b / W_Lung_b)) 
                  - (DLC_Lung * QC) * (M_Lung_b / W_Lung_b) 
                  + ((DLC_Lung * QC) * (M_Lung_t / W_Lung_t)) / P_Lung);
dxdt_M_Lung_b = R_Lung_b;
// -----------------------------------2.2 in the lung tissue compartment---------------------
double R_Lung_t = (
    (DLC_Lung * QC) * (M_Lung_b / W_Lung_b)
    - ((DLC_Lung * QC) * (M_Lung_t / W_Lung_t)) / P_Lung
    - K_up_Lung * M_Lung_t
    + K_release_Lung * M_Lung_pc
);
dxdt_M_Lung_t = R_Lung_t;

// -------------------------------------2.3 in the lung pc compartment----------------------------
double R_lung_pc = K_up_Lung * M_Lung_t - K_release_Lung * M_Lung_pc;
dxdt_M_Lung_pc = R_lung_pc;
double C_Lung = (M_Lung_b + M_Lung_t + M_Lung_pc) / W_Lung;
double C_Lung_tissue = (M_Lung_t + M_Lung_pc) / W_Lung_t;
dxdt_AUC_Lut = C_Lung_tissue*1000;
// -------------------------------------3. brain compartment--------------------------------------

// ------------------------------------3.1 brain blood compartment--------------------------------------
// PA_Brain is a factor of 10 decreased compared to other organs
double R_Brain_b = (
    Q_Brain * ((MBA / (W_Plasma * 0.2)) - (M_Brain_b / W_Brain_b))
    - (DLC_Brain * Q_Brain) * (M_Brain_b / W_Brain_b)
    + ((DLC_Brain * Q_Brain) * (M_Brain_t / W_Brain_t)) / P_Brain
);
dxdt_M_Brain_b = R_Brain_b;

// ------------------------------------3.2 brain tissue compartment--------------------------------------

double R_Brain_t = (DLC_Brain * Q_Brain) * (M_Brain_b / W_Brain_b) - ((DLC_Brain * Q_Brain) * (M_Brain_t / W_Brain_t)) / P_Brain;
dxdt_M_Brain_t = R_Brain_t;

double C_Brain = (M_Brain_b + M_Brain_t) / W_Brain;
double C_Brain_tissue = M_Brain_t/W_Brain_t;
// --------------------------------------------4. Rest of the body -----------------
// --------------------------------------------4.1 Rest of body Blood------------
double R_Rest_b = (
    Q_Rest * ((MBA / (W_Plasma * 0.2)) - (M_Rest_b / W_Rest_b)) - 
    (DLC_Rest * Q_Rest) * (M_Rest_b / W_Rest_b) + 
    ((DLC_Rest * Q_Rest) * (M_Rest_t / W_Rest_t)) / P_Rest
);
dxdt_M_Rest_b = R_Rest_b;

// --------------------------------------------4.2 Rest of body Tissue------------
double R_Rest_t = ((DLC_Rest * Q_Rest) * (M_Rest_b / W_Rest_b) - 
                    ((DLC_Rest * Q_Rest) * (M_Rest_t / W_Rest_t)) / P_Rest);
dxdt_M_Rest_t = R_Rest_t;

double C_Rest = (M_Rest_b + M_Rest_t) / W_Rest;
double C_Rest_tissue = (M_Rest_t) / W_Rest;

// ------------------------------------5. Kidney ---------------
// ------------------------------------5.1 Kidney Blood------------

double R_urine = Kurine * (M_Kidney_b / W_Kidney_b);
double R_Kidney_b = (
    Q_Kidney * ((MBA / (W_Plasma * 0.2)) - (M_Kidney_b / W_Kidney_b))
    - (DLC_Kidney * Q_Kidney) * (M_Kidney_b / W_Kidney_b)
    + ((DLC_Kidney * Q_Kidney) * (M_Kidney_t  / W_Kidney_t)) / P_Kidney
    - R_urine
);
dxdt_M_Urine = R_urine;
dxdt_M_Kidney_b = R_Kidney_b;


// ------------------------------------5.2 Kidney Tissue------------
double R_Kidney_t = (
    (DLC_Kidney * Q_Kidney) * (M_Kidney_b / W_Kidney_b)
    - ((DLC_Kidney * Q_Kidney) * (M_Kidney_t  / W_Kidney_t)) / P_Kidney
    - K_up_Kidney * M_Kidney_t
    + K_release_Kidney * M_Kidney_pc
);
dxdt_M_Kidney_t = R_Kidney_t;

// ------------------------------------5.2 Kidney pc------------
double R_Kidney_pc = K_up_Kidney * M_Kidney_t - K_release_Kidney * M_Kidney_pc;
dxdt_M_Kidney_pc = R_Kidney_pc;
double C_Kidney = (M_Kidney_b + M_Kidney_t + M_Kidney_pc) / W_Kidney;
double C_Kidney_tissue = (M_Kidney_t + M_Kidney_pc) / W_Kidney_t;
dxdt_AUC_Kt = C_Kidney_tissue*1000;
// -------------------------------------6. Spleen------
// -------------------------------------6.1 Spleen blood ------

double R_Spleen_b = (
    Q_Spleen * ((MBA / (W_Plasma * 0.2)) - (M_Spleen_b / W_Spleen_b))
    - (DLC_Spleen * Q_Spleen) * (M_Spleen_b / W_Spleen_b)
    + ((DLC_Spleen * Q_Spleen) * (M_Spleen_t / W_Spleen_t)) / P_Spleen
    - K_up_Spleen * M_Spleen_b
    + K_release_Spleen * M_Spleen_pc
);
dxdt_M_Spleen_b = R_Spleen_b;

// -------------------------------------6.2 Spleen Tissue ------
double R_Spleen_t = (
    (DLC_Spleen * Q_Spleen) * (M_Spleen_b / W_Spleen_b)
    - ((DLC_Spleen * Q_Spleen) * (M_Spleen_t / W_Spleen_t)) / P_Spleen

);
dxdt_M_Spleen_t = R_Spleen_t;

// -------------------------------------6.2 Spleen pc ------
double R_Spleen_pc = K_up_Spleen * M_Spleen_b - K_release_Spleen * M_Spleen_pc;
dxdt_M_Spleen_pc = R_Spleen_pc;
double C_Spleen = (M_Spleen_b + M_Spleen_t + M_Spleen_pc) / W_Spleen;
double C_Spleen_tissue = (M_Spleen_t + M_Spleen_pc) / W_Spleen_t;
dxdt_AUC_St = C_Spleen_tissue*1000;
// -----------------------------------7. Liver--------------------
// -----------------------------------7.1 Liver  blood --------------------


double R_Liver_b = (
    Q_Liver * (MBA / (W_Plasma * 0.2) - (M_Liver_b / W_Liver_b))
    + Q_GI * (M_GI_b / W_GI_b)
    + Q_Spleen * (M_Spleen_b / W_Spleen_b)
    - (DLC_Liver * Q_Liver) * (M_Liver_b / W_Liver_b)
    + ((DLC_Liver * Q_Liver) * (M_Liver_t / W_Liver_t)) / P_Liver
    - K_up_Liver * M_Liver_b
    + K_release_Liver * M_Liver_pc

);
dxdt_M_Liver_b = R_Liver_b;

// ---------------------------------------7.2 Liver Tissue------------
double R_bile = Kbile * (M_Liver_t / W_Liver_t);
dxdt_M_bile = R_bile;
double R_Liver_t = (
    (DLC_Liver * Q_Liver) * (M_Liver_b / W_Liver_b)
    - ((DLC_Liver * Q_Liver) * (M_Liver_t / W_Liver_t)) / P_Liver
    - R_bile
);

dxdt_M_Liver_t = R_Liver_t;


// -------------------------------------7.3 Liver pc compartment----------------
double R_Liver_pc = K_up_Liver * M_Liver_b - K_release_Liver * M_Liver_pc;
dxdt_M_Liver_pc = R_Liver_pc;
double C_Liver = (M_Liver_b + M_Liver_t + M_Liver_pc)/W_Liver;
double C_Liver_tissue = (M_Liver_t +M_Liver_pc) / W_Liver_t;

dxdt_AUC_Lt = C_Liver_tissue*1000;
// --------------------------8. GI-----------------------
// --------------------------8. 1 GI blood -----------------------
double R_GI_b = Q_GI * ((MBA / (W_Plasma * 0.2)) - (M_GI_b / W_GI_b)) 
                - (DLC_GI * Q_GI) * (M_GI_b / W_GI_b) + 
                ((DLC_GI * Q_GI) * (M_GI_t / W_GI_t)) / P_GI 
                + K_GI_b * M_GI_lumen;
                //- Kfecal * M_GI_b;

dxdt_M_GI_b = R_GI_b;

// --------------------------8. 1 GI tissue -----------------------

double R_GI_t = (
    (DLC_GI * Q_GI) * (M_GI_b / W_GI_b)
    - ((DLC_GI * Q_GI) * (M_GI_t / W_GI_t)) / P_GI
    - K_up_GI * M_GI_t
    + K_release_GI * M_GI_pc
);
dxdt_M_GI_t = R_GI_t;

// --------------------------8.2 GI pc -----------------------
double R_GI_pc = K_up_GI * M_GI_t - K_release_GI * M_GI_pc;
dxdt_M_GI_pc = R_GI_pc;
double C_GI = (M_GI_b + M_GI_t + M_GI_pc + M_GI_lumen) / W_GI;
double C_GI_tissue = (M_GI_t + M_GI_pc) / W_GI_t;
//-------------------------8.3 lumen----------
double R_GI_lumen = - K_GI_b * M_GI_lumen + R_bile - Kfecal * M_GI_lumen;
dxdt_M_GI_lumen = R_GI_lumen;
dxdt_M_fecal = Kfecal * M_GI_lumen ;

double M_Balance = (MBA + MBV + M_Lung_b + M_Lung_t + M_Lung_pc + M_Brain_b + 
                    M_Brain_t + M_Rest_b + M_Rest_t + M_Kidney_b + M_Kidney_t + 
                    M_Kidney_pc + M_Spleen_b + M_Spleen_t + M_Spleen_pc + 
                    M_Liver_b + M_Liver_t + M_Liver_pc + M_GI_b + M_GI_t + 
                    M_GI_pc + M_GI_lumen + M_bile + M_Urine +M_fecal
);

$TABLE
capture urine = M_Urine;
capture bile = M_bile;
capture feces = M_fecal;
capture Liver  = C_Liver*1000; // C_Liver with unit of mg/kg, to convert to ng/g
capture Liver_t = C_Liver_tissue*1000;
capture Kidney  = C_Kidney*1000;
capture Kidney_t = C_Kidney_tissue*1000;
capture Spleen  = C_Spleen*1000;
capture Spleen_t = C_Spleen_tissue*1000;
capture Brain  = C_Brain*1000;
capture Brain_t  = C_Brain_tissue*1000;
capture Lung  = C_Lung*1000;
capture Lung_t = C_Lung_tissue*1000;
capture GI  = C_GI*1000;
capture GI_t  = C_GI_tissue*1000;
capture Rest  = C_Rest*1000;
capture Rest_t = C_Rest_tissue*1000;
capture M_tot  = M_Balance;
capture AUC_Liver_t = AUC_Lt;
capture AUC_Kidney_t = AUC_Kt;
capture AUC_Spleen_t = AUC_St;
capture AUC_Lung_t = AUC_Lut;
capture AUC_blood = AUC_Plasma;
capture Q = Q_balance;
capture Plasma = C_Plasma*1000;
capture iv = MBV;
capture MA = MBA;
'

# -------------------------------------------------------------------------

