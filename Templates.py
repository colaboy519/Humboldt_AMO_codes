
def createHinput(r_max, k_r, N_r, max_l):
    return f'''*                                                                       *
*                  INPUT FILE FOR PROGRAM HSPH                          *
*                                                                       *
*************************************************************************
*
 BASIS SPECIFICATION
*-------------------
*
 Maximum value of r           : {round(r_max,7)}
*
 Order of the B-spline for r  : {k_r}
 Number of B-splines for r    : {N_r}
 Type of knot vector for r    : 1         ! 1 - linear
 Parameters for knot sequence : 1.05 40   ! may be two
*
 Maximum value of l           : {max_l}
*
*
 COMPUTATION SPECIFICATION
*-------------------------
*
 Type of computation         : V      ! 'N' - Compute eigenvalues only
*                                     ! 'V' - Compute eigenvalues and eigenvectors.
*
 Type of range               : A      ! 'A' - All eigenvalues will be found.
*                                     ! 'V' - All eigenvalues in the half-open interval (VL,VU] will be found.
*                                     ! 'I' - The IL-th through IU-th eigenvalues will be found.
 Value of VL                 :-1000.0 !
 Value of VU                 : 0.3    ! VL < VU
*
 Value of IL                 :  1     !
 Value of IU                 :  20    ! 1 <= IL <= IU <= N
*
*
 POTENTIAL SPECIFICATION
*-----------------------
*
 The charge of nuclear  :   1.0
*
*
***********************************************************************
'''





def bde1(R_int_,m_min_,m_max_,x_max_,N_x_,N_y_,k_x_,k_y_):
    return f'''*                                                                       *
*                  INPUT FILE FOR PROGRAM B_DAM_ECS_1e                  *
*                                                                       *
*************************************************************************
*
*
 MAIN DATA
*---------
*
 Internuclear distance       :   {round(R_int_,7)}
*
 Starting  quantum number m  :   {m_min_}
 Finishing quantum number m  :   {m_max_}
 
*
 Type     of the potential   :   1    ! See below
*
 Symmetry of the potential   :   2    ! 1: V(x,y)=V(x,-y),others: V(x,y)/=V(x,-y)
*
 State selection type        :   1    ! FOR SYMMETRIC POTENTIAL
*                                         1 - all
*                                         2 - with y-even wave functions
*                                         3 - with y-odd  wave functions
*                                         4 - with gerade symmerty
*                                         5 - with ungerade symmetry
 Use separation              :   1    !
*                                     !     1 - yes
*                                     ! other - no
 BASIS SPECIFICATION
*-------------------
*
 Maximum value of x           : {x_max_}
*
 Order of the B-spline for x  : {k_x_}
 Number of B-splines for x    : {N_x_}
 Type of knot vector for x    : 1      ! 1 - linear
 Parameters for knot sequence : 
*
 Order of the B-spline for y  : {k_y_}
 Number of B-splines for y    : {N_y_}     ! add twice for symmetric potential
 Type of knot vector for y    : 1      ! must have symmetry like potential
 Parameters for knot sequence :
*
*
 COMPUTATION SPECIFICATION
*-------------------------
*
 Type of computation         : V      ! 'N' - Compute eigenvalues only
*                                     ! 'V' - Compute eigenvalues and eigenvectors.
*
 Type of range               : A      ! 'A' - All eigenvalues will be found.
*                                     ! 'V' - All eigenvalues in the half-open interval (VL,VU] will be found.
*                                     ! 'I' - The IL-th through IU-th eigenvalues will be found.
 Value of VL                 :-1000.0 !
 Value of VU                 : 0.3    ! VL < VU
*
 Value of IL                 :  1     !
 Value of IU                 :  20    ! 1 <= IL <= IU <= N
*
*
*POTENTIAL SPECIFICATION
*-----------------------
*
 TYPE 1
*================================================================
*           Z  /   Z  /
* V(r ,r )=  a/r +  b/r  <=> U(x,y)=R( Z + Z )x + R( Z - Z )y
*    a  b    /  a   /  b                a   b         b   a
*================================================================
*
 The charge of nuclear A :   1.0
 The charge of nuclear B :   0.0
*
 TYPE 2
*================================================================
*
*================================================================
*
 !Parameter              :   1.0
 !Parameter and so on    :   1.0
*
*
*
***********************************************************************
'''


def SigmaConf(CIs_,CIp_,CId_,CIf_,m_max_):
    return f'''*                                                                            *
*         FILE WITH CONFIGURATION SERIES  FOR PROGRAM CI_Two_El              *
*                                                                            *
******************************************************************************
*
*
*
*===================================================================================
 STATES SPECIFICATION AND TYPES OF DEFINITION
*-----------------------------------------------------------------------------------
$DATA
*-----------------------------------------------------------------------------------
*                                 !
 Total angular momentum      : 0  !  
*
 Spin                        : 0  ! both - if it is neither "0" nor "1"
*                                 !
 Type of sigma states        : 2  ! Used if total angular momentum = 0.
*                                 ! 0 - plus, 1 - minus, others - both. 
*                                 !
 State symmetry type         : 2  ! FOR SYMMETRIC POTENTIAL
*                                 !       0 - with gerade   symmerty
*                                 !       1 - with ungerade symmetry
*                                 !  others - both types of symmetry
*Definition:
*
 For different spin          : 1  !Used if spin /= 0 or 1.
*                                 !Definition: 1 - the same, others: different. 
*                                 ! if 1     : ST  = ""   - for any spin value  
*                                 ! others   : ST  = "S" - for spin = 0 (singlet)
*                                 !            ST  = "T" - for spin = 1 (triplet)
*
 For different sigma states  : 1  !Used if total angular momentum = 0.
*                                 !Definition: 1 - the same, others: different. 
*                                 ! if 1     : PM = ""   - for any type  
*                                 ! others   : PM = "+"  - for plus-sigma  states
*                                 !            PM = "-"  - for minus-sigma states
*
 For different symmetry types: 0  !Used for any state symmetry type 
*                                 !Definition: 1 - the same, others: different. 
*                                 ! if 1     : GU = ""   - for any type  
*                                 ! others   : GU = "GG"  or
*                                 !            GU = "GU"  or
*                                 !            GU = "UG"  or
*                                 !            GU = "UU"  
*-----------------------------------------------------------------------------------
*
*===========================================================================
 APPROXIMATION SPECIFICATION
*---------------------------------------------------------------------------
*
 Number of l-terms      : 10      ! Exact: infinity. Normal : 15 for HeH+; 25
*for H_2
**Number of gaussian quadrature points
*
 for integral Dy        : 8       ! Exact: k_y + 1 + max_m
 for integral Lx        : 10      ! Exact: k_x + 1 + max_m
 for main integral      : 10      ! There is not exact value. Normal : 10
*
*==========================================================================
 OUTPUT SPECIFICATION
*---------------------------------------------------------------------
*   -1 ==  maximum possible
* If a number is more then maximum it used maximum.
*-----------------------------------------------------------------------
 Print energy values in log file for states with number up to : 50
*
 Print coeff. values in log file for states with number up to : 0
 Print coeff. values with absolute value more then 10^(-n), n : 0
*
 Print data in output file for states with number up to       : -1
*
*---------------------------------------------------------------------
 CONFIGURATION SERIES
*-----------------------------------------------------------------------
*
 Using restriction to quantities of nodes : 1 ! ONLY IF SEPARATION IS POSSIBLE
*                                             ! if 1   - with using restriction
*                                             ! others - without using restriction
*         
*--------------------------------------------------------------------
* With using restriction numbers below mean index numbers of states in list, 
* which contains only states satisfying conditions for quantities of nodes.
*-------------------------------------------------------------------------
* FULLCONFIG = CONFIG[M][ST PM GU]: m1 m2   m1,m2 = s,p,d,f,…     m1<= m2
* M  - Total angular momentum (<9. If it is need more, 
* it should be changed in code)
* ST - Singlet-Triplet ("S" - singlet, "T" - triplet)
* PM - Plus-Minus sigma reflection symmetry states type 
*      ("+" - plus, "-" - minus)
* GU - Gerade-Ungerade inversion symmetry states type ("GG", GU", "UG" or
*         "UU")
* EXAMPLE:| CONFIG[0][S +]: p p or | CONFIG[1][T GU]: s p or | CONFIG[1]: p d
*------------------------------------------------------------------------
*
******** S_g with l_max = {m_max_} *************
*
 CONFIG[0][GG]: s s
{CIs_}  
 CONFIG[0][UU]: s s
{CIs_}  
 CONFIG[0][GG]: p p
{CIp_}    
 CONFIG[0][UU]: p p
{CIp_}  
 CONFIG[0][GG]: d d
{CId_}    
 CONFIG[0][UU]: d d
{CId_} 
 CONFIG[0][GG]: f f
{CIf_}    
 CONFIG[0][UU]: f f
{CIf_} 
******** S_u with l_max = {m_max_} *************
*
*
 CONFIG[0][GU]: s s
{CIs_}
 CONFIG[0][UG]: s s
{CIs_}
 CONFIG[0][GU]: p p
{CIp_}  
 CONFIG[0][UG]: p p
{CIp_}
 CONFIG[0][GU]: d d
{CId_}  
 CONFIG[0][UG]: d d
{CId_}
 CONFIG[0][GU]: f f
{CIf_}  
 CONFIG[0][UG]: f f
{CIf_}
*
*
*
*
 END
************************************************************************'''


def SigmaDIP(N_DP_1_,N_DP_2_):
    return f'''*                                                                       *
*      FILE WITH INPUT DATA FOR PROGRAM DIPOLE                          *
*                                                                       *
*************************************************************************
*
*===============================================================================
$ APPROXIMATION SPECIFICATION
*-------------------------------------------------------------------------------
*
*Number of gaussian quadrature points
*
 for integral Lx        : 10       ! Exact: k_x + 1 + max_m
 for integral Qy        : 10       ! Exact: k_y + 1 + max_m
*
*=============================================================================
$ LIST OF MATRICES
*=============================================================================
*
* Example: MATRIX[S+S|s|gu] :  200  x  all
*
* 200 means use at maximum 200 states. 
* all means use all possible states
*
 MATRIX[S+S|s|gu] :  {N_DP_1_} x {N_DP_2_}
*
*
$ END
*
***********************************************************************'''


def LaserPulse(pulseenergy_,pulseduration_,laseramplitude_):
    return f'''*                                                                       *
*                  INPUT FILE DESCRIBING PULSE PARAMETERS               *
*                                                                       *
*************************************************************************
*
$ ENVELOPE BASIS
*-------------------
*
 Physical quantity : A   ! "F" (electric field) or "A" (vector potential)
*
$ CARRIER-ENVELOPE FREQUENCY
*----------------------------
*
 Type :   F    ! "V" - varying wavelength / frequency ( numerical value is ignored)
*              ! "F" - fixed wavelength / frequency   ( numerical value is used) 
*              ! "S" - static (next records are ignored)
*
 Units:   w    ! "L" - the wavelength   is specified in nm
*              ! "l" - the wavelength   is specified in a.u.
*              ! "W" - the frequency    is specified in eV
*              ! "w" - the frequency    is specified in a.u.
*              ! "T" - the cycle period is specified in fs
*              ! "t" - the cycle period is specified in a.u.
*
 Numerical value  :  {pulseenergy_}
*
$ SIDES SPECIFICATION
*------------------- 
*
 Shape :  C    !   "G" : Gaussian pulse envelope 
*              !   "S" : Sech pulse envelope 
*              !   "C" : Cos2 pulse envelope 
*              !   "P" : Parabolic pulse envelope 
*              !   "L" : Linear ascent and descent. 
*              !   "N" : no sides (pure flat pulse, next records are ignored)
*
 Type of duration definition:  s   ! "u" - united definition
*                                  ! "s" - specific definition
*
 Type :   F    ! "V" - varying sides duration ( numerical value is ignored)
*              ! "F" - fixed sides duration   ( numerical value is used)
*
 Units :  t    ! "C" - the side duration is specified in numbers of cycles
*              ! "T" - the side duration is specified in fs
*              ! "t" - the side duration is specified in a.u.
*
 Numerical value  : {pulseduration_}
*
 Parameter :                 ! If it is not empty, then it substitutes
*                            ! the default value of the parameter 
*                            ! restricting the duration of an infinite pulse
$ FLAT DURATION
*------------------- 
*
 Type :   N    ! "V" - varying flat duration ( numerical value is ignored)
*              ! "F" - fixed flat duration   ( numerical value is used) 
*              ! "N" - no flat (next records are ignored)
*
 Units:   C    ! "C" - the flat duration is specified in the numbers of cycles
*              ! "T" - the flat duration is specified in fs
*              ! "t" - the flat duration is specified in a.u.
*
 Numerical value  :      
*
*
$ PEAK CHARACTERISTIC
*------------------- 
*
 Type :  F     ! "V" - varying peak characteristic ( numerical value is ignored)
*              ! "F" - fixed peak characteristic   ( numerical value is used)
*
 Units : I     ! "I" - the peak intensity        is specified in W/cm^2
*              ! "F" - the peak electric field   is specified in a.u.
*              ! "A" - the peak vector potential is specified in a.u.
*
 Numerical value  : {laseramplitude_[0]}D+{laseramplitude_[1]}
*
$ CARRIER-ENVELOPE PHASE
*----------------------- 
*
* phi = x * pi/8
*
 x :   4       !
*
$ POLARISATION
*-----------------------
*
 Ellipticity : 0.0   !
*'''


def Tp2input(input1e_,input2e_,inputdip_,NumberofSymmetries_):
    return f'''*                                                                       *
*                  INPUT FILE FOR PROGRAM TP_BDE_2e                     *
*                                                                       *
*************************************************************************
*
 BASIS SPECIFICATION
*-------------------
*
 The name of the system             : H_H
 The name of input data for 1e      : {input1e_}
 The name of input data for 2e      : {input2e_}
 The name of input data for dipole  : {inputdip_}
 The number of used symmetries      : {NumberofSymmetries_}
*
*-----------------------------------------------------------------
* Specification of names for dipole matrices between symmetries.
* Without the extension. E.g. : S_P_s_gu
*-----------------------------------------------------------------
*
 Symmetry1 - Symmetry2  : SpS_s_gu
* Symmetry2 - Symmetry3  : P_D_s_ug
* Symmetry3 - Symmetry4  : D_F_s_gu
* Symmetry4 - Symmetry5  : F_G_s_ug
* Symmetry5 - Symmetry6  : G_H_s_gu
* Symmetry6 - Symmetry7  : H_I_s_ug
* Symmetry7 - Symmetry8  : I_J_s_gu
*
*
***********************************************************************'''

def nucfixbase(Jmin_,Jmax_,Rmax_,Bk_,Nk_,r_red_):
    return f'''*
*
*        Input file for calculating the nuclear 
*        wave functions using a B-spline basis:
*
*
*    Data concerning the molecular system:
*    -------------------------------------
*
   Masses m_A and m_B (or mu_(reduced$) and -1.0D+00)   :   {r_red_}D+00  -1.0D+00
*
   Range of rotational quantum numbers J (begin, end)  :   {Jmin_}  {Jmax_}
*
   Electronic angular and spin momenta (projected on z):   0  0
*
*
*    Data concerning the B-spline basis:
*    -----------------------------------
*
   R_(max) (=box radius)    :    {Rmax_}
*
   Order of the B-splines   :    {Bk_}
*
   Number of B-splines      :    {Nk_}
*
   Type of knot sequence    :    0      !  0: linear, 1: sine-like
*
   First non-zero knot point: 1.0D-8   ! (Not used for linear knot sequence)  
*
*
*
*    Data concerning the number of eigenvalues and -vectors to be stored:
*    --------------------------------------------------------------------
*
   Save parameter                           :   2   ! 0: Do not store
*                                                   !    anything.
*                                                   ! 1: Store eigenvalues.
*                                                   ! 2: Store eigenvalues
*                                                   !    and -vectors.
*
   Range of states to be stored (begin, end):   1  5000    ! These parameters
*                                                       ! are ignored, if
*                                                       ! the save parameter
*                                                       ! is set to 0.'''