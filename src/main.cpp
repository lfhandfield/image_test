#include "primitive.h"
#define REGRESS_POWER (8)


using namespace LFHPrimitive;


void myfilestreamcopy(FILE *out, FILE *in, unsigned int size, char* buffer, int buffer_magnitude){
	int bsize = 1 << buffer_magnitude;
	while(size >> buffer_magnitude){
		fread(buffer, sizeof(char),bsize,in);
		fwrite(buffer, sizeof(char),bsize,out);
		size -= bsize;
	}
	fread(buffer, sizeof(char),size,in);
	fwrite(buffer, sizeof(char),size,out);
}

void myfilecopy(FILE *out, FILE *in, char* buffer, int buffer_magnitude){
	int bsize = 1 << buffer_magnitude;
	int sout = bsize;
	while(sout == bsize){
		sout = fread(buffer, sizeof(char),bsize,in);
		if (sout > 0) fwrite(buffer, sizeof(char),sout,out);
	}
}

double cleveratof(char* buf){
	double tmp = atof(buf);
	return(tmp != 0.0f ? tmp : (buf[0] != '0' ? NAN : 0.0f));
}

void renorm(double* arr, int arrsize){
	double sum=0.0f;
	int i;
	int l=arrsize;
	for(i=0;i<arrsize;i++){
		if (isnan(arr[i])||isinf(arr[i])) l--;
		else sum += arr[i];
	}
	if (l > 0) {
		sum /= l;
		for(i=0;i<arrsize;i++) arr[i] -= sum;
	}
}


double tuple_distance(const Tuple<double, 3> &a,const Tuple<double, 3> &b ){
	double tmp = a[0] -b[0];
	double tmp2 = a[1] -b[1];
	return(tmp*tmp + tmp2*tmp2);
	
	}

double tuple_weight(const Tuple<double, 3> &a ){
	return(a[2]);	
	}

	
	 
	
	
	/*    
	    
	  
	  
	 
	
	[2 a f_1f_3b f_1f_4b f_1f_5b] [F_1*X_1 + M]    [F_1*Y_1 + K]  
	[a 2 f_2f_3b f_2f_4b f_2f_5b] [F_2*X_2 + M]    [F_2*Y_2 + K]  
	[f_1f_3b f_2f_3b 2 c d]       [F_3*X_3 + N] =  [F_3*Y_3 + L]  
	[f_1f_4b f_2f_4b c 2 e]       [F_4*X_4 + N]    [F_4*Y_4 + L]  
	[f_1f_5b f_2f_5b d e 2]       [F_5*X_5 + N]    [F_5*Y_5 + L]  
	
	[2 a] [X_1] = [Y_1]
	[a 2] [X_2]   [Y_2]

	[2 c d] [X_3]   [Y_3]
	[c 2 e] [X_4] = [Y_4]
	[d e 2] [X_5]   [Y_5]
	
	
	
	
	2X_1 + aX_2 = Y_1 
	aX_1 + 2X_2 = Y_2 

	2X_3 + cX_4 + dX_5 = Y_3 
	cX_3 + 2X_4 + eX_5 = Y_4 
	dX_3 + eX_4 + 2X_5 = Y_5 
	
	2*N * (3 + c + d + e) + b*(f_1 + f_2)*M*(f_3 + f_4 + f_5)    
	F_3*X_3*(2+c+d) + F_4*X_4*(2+c+e) + F_5*X_5*(2+e+d) + b*(F_1*X_1f_1 + F_2*X_2f_2)*(f_3 + f_4 + f_5)     =  (F_3*Y_3 + F_4*Y_4 + F_5*Y_5) + 3L
	

	2*N * (3 + c + d + e) + b*(f_1 + f_2)*M*(f_3 + f_4 + f_5)    
	b*(f_3 + f_4 + f_5)*(F_1*X_1f_1 + F_2*X_2f_2)   = 3K

	2*M * (2 + a) + b*N*(f_1 + f_2)*(f_3 + f_4 + f_5)    
	b*(F_3*X_3f_3 + F_4*X_4f_4 + F_5*X_5f_5)*(f_1 + f_2)   = 2L
	
	[2*(2 * a)              b*(f_1 + f_2)*(f_3 + f_4 + f_5)] [M] = [2L]
	[b*(f_1 + f_2)*(f_3 + f_4 + f_5)      2*(3 + c + d + e)] [N]   [3K]
	
	
	F_1*Y_1 + K
	F_2*Y_2 + K
	F_3*Y_3 + L
	F_4*Y_4 + L
	F_5*Y_5 + L
	
	2*(2 * a)*M + b*(f_1 + f_2)*(f_3 + f_4 + f_5)*N          = 3*L  = 2*(F_1*Y_1*f_1  + F_2*Y_2*f_2)/(1 - f_1 - f_2)
	b*(f_1 + f_2)*(f_3 + f_4 + f_5)*M + 2*(3 + c + d + e)*N  = 3*K  = 3*(F_3*Y_3*f_3  + F_4*Y_4*f_4+ F_5*Y_5*f_5)/(1 - f_3 - f_4 - f_5)
	
	b*(f_3 + f_4 + f_5)*(F_1*X_1f_1 + F_2*X_2f_2)   = 3K

		
	F_1*X_1f_1 +  = -F_2*X_2f_2
	F_3*X_3f_3 + F_4*X_4f_4 + F_5*X_5f_5 = 0
	
	2*N * (3 + c + d + e) + b*(f_1 + f_2)*M*(f_3 + f_4 + f_5)    
	F_3*X_3*(2+c+d) + F_4*X_4*(2+c+e) + F_5*X_5*(2+e+d) + b*(F_1*X_1f_1 + F_2*X_2f_2)*(f_3 + f_4 + f_5)     =  (F_3*Y_3 + F_4*Y_4 + F_5*Y_5) + 3L

		
	suppose F_i = 1 / f_i:
	X_1 + X_2 =0
	X_3 + X_4 + X_5 =0

	
	Y_1 = -Y_2 

	2X_3 + cX_4 + d(-X_3 -X_4) = Y_3 
	cX_3 + 2X_4 + e(-X_3 -X_4) = Y_4 
	dX_3 + eX_4 + 2(-X_3 -X_4) = Y_5 

	then!
	Y_3 + Y_4 + Y_5 = 0
	
	
	
	SO!
	
	Given Z_1 Z_2 Z_3 Z_4 Z_5
	
	L = (f_1 *Z_1 + f_2*Z_2) / (f_1 + f_2)
	F_i = 1 / (f_i); 
	implies Y_1 + Y_2 = 0
	
	F_1Y_1 + L = Z_1
	F_1Y_1*f_1 + L*f_1 = Z_1*f_1
	F_2Y_2*f_2 + L*f_2 = Z_2*f_2
	
	Y_1+ Y_2 + L*(f_1 +f_2) = Z_1*f_1 + Z_2*f_2  
	Y_1+Y_2 = 0
	
	which then implies:
	
	X_1 + X_2 = 0
	
	
	which then implies that the linear system of equation is solved!
	
	
	
	Y_1 = f_1*(Z_1 - L) 
	
     
	
	Basically: 
		
	[2 a f_1f_3b f_1f_4b f_1f_5b] [M + X_1/f_1]    [K + Y_1/f_1]  
	[a 2 f_2f_3b f_2f_4b f_2f_5b] [M + X_2/f_2]    [K + Y_2/f_2]  
	[f_1f_3b f_2f_3b 2 c d]       [N + X_3/f_3] =  [L + Y_3/f_3]     
	[f_1f_4b f_2f_4b c 2 e]       [N + X_4/f_4]    [L + Y_4/f_4]  
	[f_1f_5b f_2f_5b d e 2]       [N + X_5/f_5]    [L + Y_5/f_5]  

	[2 a] [X_1] = [Y_1]
	[a 2] [X_2]   [Y_2]

     
	
	[2 c d] [X_3]   [Y_3] 
	[c 2 e] [X_4] = [Y_4]
	[d e 2] [X_5]   [Y_5]

	[2*(2 * a)              b*(f_1 + f_2)*(f_3 + f_4 + f_5)] [M] = [2L]
	[b*(f_1 + f_2)*(f_3 + f_4 + f_5)      2*(3 + c + d + e)] [N]   [3K]
	
	 
	example time!
	
	    
	[4 2 1]   =   [2]
	[2 4 1] X =   [5]
	[1 1 4]   =   [9]
	
	
	
			
	[2 a f_1f_3b f_1f_4b f_1f_5b] [M + X_1/f_1]    [K + Y_1/f_1]  
	[a 2 f_2f_3b f_2f_4b f_2f_5b] [M + X_2/f_2]    [K + Y_2/f_2]  
	[f_1f_3b f_2f_3b 2 c d]       [N + X_3/f_3] =  [L + Y_3/f_3]     
	[f_1f_4b f_2f_4b c 2 e]       [N + X_4/f_4]    [L + Y_4/f_4]  
	[f_1f_5b f_2f_5b d e 2]       [N + X_5/f_5]    [L + Y_5/f_5]  
	
	(4+2a)M  +N *... (2+a) * X_1/f_1 (2+a) * X_2/f_1 + (f_1+f_2)b*(X_1 + X_2 + X_3) =  2K + Y_1/f_1 + Y_2/f_2
	
	 X_1/f_1 (2+a) * X_2/f_1 + (f_1+f_2)b*(X_1 + X_2 + X_3) =   Y_1/f_1 + Y_2/f_2
	
	  
	*/
 
		
class bobo{
	public:     
	
	static void addsix(int &a, const int &b){a = b + 6;}
	
	};
	
class stupid{
	public:
	int haha;
	stupid(){}
	stupid(int kk):haha(kk){}
	void memmove(stupid& source){
		haha = source.haha;
		source.haha =0;
		printf("moved %i, from %X to %X\n",haha, &source, this);
		}
	void show(FILE *f){printf("%i", haha);}
	bool operator<(const stupid& ot)const{return haha < ot.haha;}
	bool operator<=(const stupid& ot)const{return haha <= ot.haha;}
	bool operator>(const stupid& ot)const{return haha > ot.haha;}
	bool operator>=(const stupid& ot)const{return haha >= ot.haha;}
	bool operator==(const stupid& ot)const{return haha == ot.haha;}
	bool operator!=(const stupid& ot)const{return haha != ot.haha;}
	
	typedef YESNO< false > IsPOD;
	void show(FILE*f , int level)const {fprintf(f,"%i", haha);}
	
	};
	
int main (int argc, char * const argv[]) {
	
	DataGrid<double ,  2> fun;
	
	unsigned int dim[2];
	
	dim[0] = 10;
	dim[1] = 10;
	fun.setSizes(dim);
	
	for(dim[0]=0;dim[0]<10;dim[0]++)
	for(dim[1]=0;dim[1]<10;dim[1]++)
	fun(dim) = (-5 + ((int)dim[1])) * ((double)dim[0]);
	
	
	DataGrid<double , 1> slice = fun(7);
	  
	ExOp::show(fun);
	ExOp::show(slice);
	
	printf("%f\n", slice(1));
	printf("%f\n", fun(1,8));
	
	
	//slice(5) =1.56f;
	ExOp::show(slice);
	 
	exit(0);
	 
	
	RBTofDoom<stupid> brbr,trtr;
	
for(int ff = 0; ff< 7;ff++) brbr.insert(stupid(rand()));

brbr.show();
	brbr.remove(stupid(41));

 	 brbr.show();
	 
	 brbr.remove(stupid(634));

	 brbr.show();        
	 brbr.remove(stupid(6334));

	 brbr.show();
	  	 brbr.remove(stupid(18467));

	 brbr.show();
	 
	 	    brbr.remove(stupid(26500));

	 brbr.show();
	 
	 
	exit(0);	  
	
	FILE* f = fopen("./tmptmp", "wb+");
	brbr.save(f);
	fclose(f);
	
	f = fopen("./tmptmp", "rb+");
	trtr.load(f);
	fclose(f);
	
	trtr.show();
	
	exit(0);  
	
	
	/*	
	
	SpacePartition<unsigned int, 3,0, int> spacep;
	
	Tuple<unsigned int, 3> min[2];   
	min[0][0] = 0; min[0][1] = 0; min[0][2] = 0;   
	min[1][0] = 0; min[1][1] = 15; min[1][2] = 0;     
	
	
	spacep.insert(HyperPosition<unsigned int, 3,0>( min[0], 4 , 2),0);
	min[0][0] = 0; min[0][1] = 15; min[0][2] = 0;     
	spacep.insert(HyperPosition<unsigned int, 3,0>( min[0], 4 , 2),0);
 	min[0][0] = 0; min[0][1] = 0; min[0][2] = 15;     
	spacep.insert(HyperPosition<unsigned int, 3,0>( min[0], 4 , 2),0);
	min[0][0] = 0; min[0][1] = 15; min[0][2] = 15;     
	spacep.insert(HyperPosition<unsigned int, 3,0>( min[0], 4 , 2),0);
	    
	min[0][0] = 256; min[0][1] = 15; min[0][2] = 15;      
	HyperPosition<unsigned int, 3,0> funfun = HyperPosition<unsigned int, 3,0>( min[0], 4 , 2);
	int i;
	
//	Tuple <double, 3>  center; 
	 
//	vector<int> intersection
//	spacep.hypersphereIntersection(center,2048.0f, intersection);
	
 
	
	
	//spacep
	
	
*/	
	
//	unsigned int coor[] = { 15, 15,15,8,6,4};
//	HyperBox tmpbox = HyperBox::makeFromCorner( coor,  coor+3);
	  
	/*
	   
	Vector<Tuple<double,3 > > coors; 

  		 
	DataGrid<double, 2, LFHVECTOR_NORMAL> predata;
	DataGrid<double, 2, LFHVECTOR_NORMAL> postdata;
	   
	unsigned int dims[2];
	unsigned int coor[2];
	Tuple<double,3 > icoor;
	 
	dims[1] = 40; // nbpts              
	            
	dims[0] = 3; 
	       
	predata.setSizes(dims);     
	postdata.setSizes(dims);     
	 
	              
	for(coor[1] =0;coor[1] <dims[1];coor[1]++){
//		icoor[0] = ((double)rand()) *(1.0f / RAND_MAX);
//		icoor[1] = ((double)rand()) *(1.0f / RAND_MAX);  


		icoor[0] = ((double)(rand() % 100)) *0.01f;
		icoor[1] = ((double)(rand() % 100)) *0.01f;  
//		icoor[2] = ((double)(rand() % 100)) *0.01f;  
		icoor[2] =  ((double)(rand() % 100)) *0.01f;     
		icoor[2] = 1.0f;
		coors.push_back(icoor);      
		for(coor[0] =0;coor[0] <dims[0];coor[0]++){    
		predata(coor) = ((double)rand()) *(255.0f / RAND_MAX);  
		postdata(coor) = 0.0f;  
		}
		
		}
			
	//   
	LFHPrimitive::GaussianProcess::GPweightsregression(postdata, predata, coors , tuple_distance, tuple_weight);

	
	*/
	
	
	
	
	
	return 0;
}
