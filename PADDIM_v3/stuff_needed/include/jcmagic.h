/** for the 2D stuff **/
#define JC2D   76

#define JCGRID 10    /** 2D Grid **/
#define JCTEMP 11    /** 2D Temp **/
#define JCCHEM 12    /** 2D Chem **/
#define JCSTREAM 13  /** 2D stream function : psi, d psi/dx , - d psi/ dy **/
#define JCPSI 14     /** 2D stream function : psi, d psi/dx , - d psi/ dy **/
#define JCPSIDX 15   /** 2D stream function : psi, d psi/dx , - d psi/ dy **/
#define JCPSIDY 16   /** 2D stream function : psi, d psi/dx , - d psi/ dy **/


#define JC3D   77
#define JC3R   78
#define JC3S   79

#define JPTEMP 20   /*** 3D temp field **/
#define JPAX 21     /*** 3D stream x field **/
#define JPAY 22     /*** 3D stream y field **/
#define JPVELX 23   /*** 3D  vel x field **/
#define JPVELY 24   /*** 3D vel y  field **/
#define JPVELZ 25   /*** 3D vel z  field **/
#define JPPRES 26   /*** 3D pressure field **/
#define JPCHEM 27   /*** 3D chem field **/
#define JPBX 28   /*** 3D BX field **/
#define JPBY 29   /*** 3D BY field **/
#define JPBZ 30   /*** 3D BZ field **/

#ifdef JCBYNAME
static char jcnamebytype[40][8] = {"","","","","","","","","","",
				"Grid","Temp","Chem","Stream","PSI","PSI/DX","PSI/DY","","","",
				"Temp","AX","AY","VX","VY","VZ","P","C","BX","BY",
				"BZ","","","","","","","","","",
				};
#endif



#define JPGRID 37   /*** 3D grid field **/

#ifdef NOFORT
#ifdef __cplusplus
extern "C" {
#endif
int jcopen(int ,int ,char *,char *);
int jcwinfo(int ,float *,float *,int *,float *,int *,int *, int *);
int jcwrite(int , int , int *,float *,float *,float *,int *,int *,int *,int * );
int jcskip(int , int , int *,float *);
int jcrinfo(int ,float *,float *,int *,float *,int *,int *, int *);
int jcread(int ,int *,float *,float *,float *,int *,int *,int *);
int jcclose(int );

#ifdef __cplusplus
}
#endif

#endif




