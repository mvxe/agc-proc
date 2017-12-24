/*
    alpha-gamma-time processing
    Copyright (C) 2017  Mario Vretenar

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <cstdio>
#include <cstdlib>
#include <string>
#include <fstream>
#include <cmath>

using namespace std;

int alpha_zero,gamma_zero;      	//the zero value
int alpha_thresh,gamma_thresh;  	//the threshold value
unsigned alpha_Nmes,gamma_Nmes; 	//the total number of measurements per file (alpha.dat and gamma.dat)
unsigned time_len;
unsigned alpha_Nch,gamma_Nch;		//the total number of measurements (time.dat) will be [alpha_Nch]x[gamma_Nch]x[time_len]
unsigned step_alpha,step_gamma;
unsigned alpha_maxN, gamma_maxN;
unsigned alpha_max=0;
unsigned gamma_max=0;
		
void _help()
{
	printf ("Arguments: <path> <mode>\n"
	        "<path>: path to folder containing measurements\n"
	        "<mode>: 0 - integration on time resolved\n"
	        "        1 - integration on time resolved minus background\n"
	        "        2 - peak - background ratio\n"
	        "        3 - integration on time resolved minus background (background calculated using total num of events)\n"
	        "        4 - peak - background ratio (background calculated using total num of events)\n");
}

double xmax,xmin,ymax,ymin,xstep,ystep;
void _load_conf(string fpath)
{
	ifstream t((fpath + "/agc_conf.txt").c_str());
	string conffile((istreambuf_iterator<char>(t)), istreambuf_iterator<char>());            //put the conf file into a string
	if (conffile.size()==0) {printf("Missing agc_conf.txt. Aborting...\n");exit(0);}
	
	char tmpch; int z;
	bool alpha_edge,gamma_edge;
	size_t pos_alpha_thresh = conffile.find("alpha_thresh(-8192 - 8191):");
		if (pos_alpha_thresh != string::npos){
			pos_alpha_thresh+=27;
			sscanf(conffile.substr(pos_alpha_thresh).c_str(), "%d", &alpha_thresh);
			printf("alpha_thresh=%d\n",alpha_thresh);
		}else {printf("Error in alpha_thresh.\n"); exit(0);}
	size_t pos_alpha_edge = conffile.find("alpha_edge(Rising (R) or Falling (F)):");
		if (pos_alpha_edge != string::npos){
			pos_alpha_edge+=38;
			z=0; do {sscanf(conffile.substr(pos_alpha_edge+z).c_str(), "%c", &tmpch); z++;}
			while (isspace(tmpch));
			if (tmpch=='R') alpha_edge=0;
			else if (tmpch=='F') alpha_edge=1;
			else {printf("Error in alpha_edge. Must be F or R!\n"); exit(0);}
		}else {printf("Error in alpha_edge.\n"); exit(0);}
	size_t pos_gamma_thresh = conffile.find("gamma_thresh(-8192 - 8191):");
		if (pos_gamma_thresh != string::npos){
			pos_gamma_thresh+=27;
			sscanf(conffile.substr(pos_gamma_thresh).c_str(), "%d", &gamma_thresh);
			printf("gamma_thresh=%d\n",gamma_thresh);
		}else {printf("Error in gamma_thresh.\n"); exit(0);}
	size_t pos_gamma_edge = conffile.find("gamma_edge(Rising (R) or Falling (F)):");
		if (pos_gamma_edge != string::npos){
			pos_gamma_edge+=38;
			z=0; do {sscanf(conffile.substr(pos_gamma_edge+z).c_str(), "%c", &tmpch); z++;}
			while (isspace(tmpch));
			if (tmpch=='R') gamma_edge=0;
			else if (tmpch=='F') gamma_edge=1;
			else {printf("Error in gamma_edge. Must be F or R!\n"); exit(0);}
		}else {printf("Error in gamma_edge.\n"); exit(0);}
	size_t pos_alpha_zero = conffile.find("alpha zero level (not needed by program, for reference):");
		if (pos_alpha_zero != string::npos){
			pos_alpha_zero+=56;
			sscanf(conffile.substr(pos_alpha_zero).c_str(), "%d", &alpha_zero);
			printf("alpha_zero=%d\n",alpha_zero);
		}else {printf("Error in alpha_zero.\n"); exit(0);}
	size_t pos_gamma_zero = conffile.find("gamma zero level (not needed by program, for reference):");
		if (pos_gamma_zero != string::npos){
			pos_gamma_zero+=56;
			sscanf(conffile.substr(pos_gamma_zero).c_str(), "%d", &gamma_zero);
			printf("gamma_zero=%d\n",gamma_zero);
		}else {printf("Error in gamma_zero.\n"); exit(0);}
	size_t pos_step_alpha = conffile.find("Time resolved alpha amplitude step:");
		if (pos_step_alpha != string::npos){
			pos_step_alpha+=35;
			sscanf(conffile.substr(pos_step_alpha).c_str(), "%u", &step_alpha);
		}else {printf("Error in step_alpha.\n"); exit(0);}
	size_t pos_step_gamma = conffile.find("Time resolved gamma amplitude step:");
		if (pos_step_gamma != string::npos){
			pos_step_gamma+=35;
			sscanf(conffile.substr(pos_step_gamma).c_str(), "%u", &step_gamma);
		}else {printf("Error in step_gamma.\n"); exit(0);}
	double interval;
	size_t pos_interval = conffile.find("trigger event(0 - 34.3597)(in seconds):");
		if (pos_interval != string::npos){
			pos_interval+=39;
			sscanf(conffile.substr(pos_interval).c_str(), "%lf", &interval);
		}else {printf("Error in interval. Delete file to regenerate from template.\n"); exit(0);}
	time_len=(unsigned)(interval*125000000);
	size_t pos_alpha_max = conffile.find("alpha_max(R edge: alpha_thresh < x < 8191, F edge: -8192 < x < alpha_thresh):");
		if (pos_alpha_max != string::npos){
			pos_alpha_max+=77;
			sscanf(conffile.substr(pos_alpha_max).c_str(), "%d", &alpha_maxN);
		}
		else if(!alpha_edge)	//we dont error on this for backwards compatibility
			alpha_maxN = 8191;
		else alpha_maxN = -8192;
	size_t pos_gamma_max = conffile.find("gamma_max(R edge: gamma_thresh < x < 8191, F edge: -8192 < x < gamma_thresh):");
		if (pos_gamma_max != string::npos){
			pos_gamma_max+=77;
			sscanf(conffile.substr(pos_gamma_max).c_str(), "%d", &gamma_maxN);
		}
		else if(!gamma_edge)	//we dont error on this for backwards compatibility
			gamma_maxN = 8191;
		else gamma_maxN = -8192;
	
	if (!alpha_edge) alpha_Nmes = alpha_maxN-alpha_thresh+1;
	else alpha_Nmes = -(alpha_maxN-alpha_thresh)+1;
	alpha_Nch = alpha_Nmes/step_alpha+1;
	if (!gamma_edge) gamma_Nmes = gamma_maxN-gamma_thresh+1;
	else gamma_Nmes = -(gamma_maxN-gamma_thresh)+1;
	gamma_Nch = gamma_Nmes/step_gamma+1;
	printf("alpha_Nmes=%u\n",alpha_Nmes);
	printf("gamma_Nmes=%u\n",gamma_Nmes);
	printf("alpha_Nch=%u\n",alpha_Nch);
	printf("gamma_Nch=%u\n",gamma_Nch);
	printf("time_len=%u\n",time_len);
	
	if (conffile.find("Observed interval before and after trigger event(0 - 34.3597)(in seconds):") != string::npos) time_len*=2;	//backward(forward) compatibility
	
	xmin=0;
	xmax=alpha_Nch-1;
	ymin=0;
	ymax=gamma_Nch-1;
	xstep=1000;
	ystep=1000;
}

FILE *Fgplot;
void _gnuplot_set(int stage, int mode)
{
	double amax;
	if(mode>2) mode-=2;
	int N=0;
	switch (stage){
	case 0: Fgplot = popen("gnuplot -persistent", "w"); //we want to see error messages
		fprintf ( Fgplot,	"set terminal wxt size 400,400\n"
					"set multiplot\n");
	    	break;
	case 1: amax=(double)alpha_max/4;
		while(amax>10) {amax/=10;N++;}
		amax=floor(amax);
		fprintf ( Fgplot,	"set tics font \"Helvetica,10\"\n"
					"set xlabel font \"Helvetica,13\"\n"
					"set ylabel font \"Helvetica,13\"\n"
					"set zlabel font \"Helvetica,13\"\n"
					"set key font \"Helvetica,10\"\n"
					"set size 0.74*0.91,0.21\n"
					"set bmargin 0\n"
					"set lmargin 0\n"
					"set rmargin 0\n"
					"set tmargin 0\n"
					"set origin 0.26+0.08*0.74,0.78\n"
					"set xrange [%lf:%lf]\n"
					"set yrange [0:%lf]\n"
					"set xtics in\n"
					"set ytics in\n"
					"set xtics format \"%%g\"\n"
					"set ytics format \"%%g\"\n"
					"set xtics offset 0,0.55\n"
					"set xtics %lf\n"
					"set xtics scale 1\n"
					"set ytics scale 1\n"
					"set ytics (%lf,%lf,%lf,%lf,%lf)\n"
					"set label '0' at screen 0.26+0.08*0.74-0.04,screen 0.78 font \"Helvetica,10\"\n"	//add a zero tic
					"set format y \"%%.1t*10^{%%T}\"\n"
					"unset xlabel\n"
					"set key left top\n"
					"set style fill transparent solid 0.5\n",(xmin+0.08*(xmax-xmin))*step_alpha,xmax*step_alpha,
						5.5*amax*pow(10,N),xstep,amax*pow(10,N),2*amax*pow(10,N),3*amax*pow(10,N),4*amax*pow(10,N),5*amax*pow(10,N));
	    	break;
	case 2: amax=(double)gamma_max/4;
		while(amax>10) {amax/=10;N++;}
		amax=floor(amax);
		fprintf ( Fgplot,	"set size 0.21,0.91*0.74\n"
					"set bmargin 0\n"
					"set lmargin 0\n"
					"set rmargin 0\n"
					"set tmargin 0\n"
					"set origin 0.01,0.01\n"
					"set yrange [%lf:%lf]\n"
					"set xrange [%lf:0] reverse\n"
					"set xtics format \" \"\n"
					"set ytics format \"%%g\"\n"
					"set ytics %lf\n"
					"set ytics rotate\n"
					"set ytics offset 14.6,0\n"
					"set xtics (%lf,%lf,%lf,%lf,%lf)\n"
					"set x2tics (%lf,%lf,%lf,%lf,%lf)\n"
					"set label '0' at screen 0.22,screen 0.92*0.74+0.02 rotate font \"Helvetica,10\"\n"	//add a zero tic
					"set format x2 \"%%.1t*10^{%%T}\"\n"
					"set x2label \"counts\\nper\\nchannel\"\n"
					"set x2label offset -1,0\n"
					"set x2tics rotate\n"
					"unset xlabel\n"
					"unset ylabel\n"
					"set label 'gamma' at screen 0.05,screen 0.05 rotate font \"Helvetica,10\"\n"		//gotta do this since gnuplot doesnt support key rotation
					"set style rectangle back fill border lc rgb \"#74c22b\"\n"	//it wont let me do transparent and border at the same time for some reason
					"set object rectangle from screen 0.038,screen 0.315-0.145 to screen 0.064,screen 0.395-0.145\n"
					"set style rectangle back fill transparent solid 0.5 fc rgb \"#74c22b\"\n"
					"set object rectangle from screen 0.038,screen 0.315-0.145 to screen 0.064,screen 0.395-0.145\n"
					"set style fill transparent solid 0.5\n",ymin*step_gamma,ymax*step_gamma-0.08*(ymax-ymin)*step_gamma,
						5.5*amax*pow(10,N),ystep,amax*pow(10,N),2*amax*pow(10,N),3*amax*pow(10,N),4*amax*pow(10,N),5*amax*pow(10,N),
						amax*pow(10,N),2*amax*pow(10,N),3*amax*pow(10,N),4*amax*pow(10,N),5*amax*pow(10,N));
	    	break;
	case -2: fprintf ( Fgplot,	"unset label\n"
					"set size 1,1\n"
					"set bmargin screen 0.01\n"
					"set lmargin screen 0.26\n"
					"set rmargin screen 0.99\n"
					"set tmargin screen 0.74\n"
					"set origin 0,0\n"
					"unset x2label\n"
					"unset x2tics\n"
					"set xrange [%lf:%lf]\n"
					"set yrange [%lf:%lf]\n"
					"set cbrange [%c:*]\n"
					//"set logscale cb\n"
					"unset xlabel\n"
					"unset ylabel\n"
					"unset zlabel\n"
					"set xtics out\n"
					"set ytics out\n"
					"set xtics scale 0.3\n"
					"set ytics scale 0.3\n"
					"unset ztics\n"
					"set xtics format \" \"\n"
					"set ytics format \" \"\n"
					"set xtics %d\n"
					"set ytics %d\n"
					"set ytics nomirror\n"
					"set pm3d map\n"
					"unset colorbox\n"
	    				"unset dgrid3d\n"
	    				"unset style fill\n"
	    				"set style rectangle back fill border lc rgb \"black\"\n"
					"set object rectangle from screen 0.26,screen 0.01 to screen 0.99,screen 0.74 fill transparent solid 1 fc rgb \"black\" behind\n"
					"set object rect from screen 0.26,screen 0 to screen 1,screen 0.01 fc rgb \"white\" fs solid 1 noborder\n"
					,xmin,xmax,ymin,ymax,(mode==2)?'*':'0',(int)(xstep/step_alpha),(int)(ystep/step_gamma));    
	    		if (mode==0 || mode==1) 
	    			fprintf ( Fgplot,"A=1.7\n"
					"set palette model XYZ functions gray**(A*0.35), gray**(A*0.5), gray**(A*0.8)\n");
			else if (mode==2)
				fprintf ( Fgplot,"set palette maxcolors 100\n"
					"set palette defined (0 \"blue\", 99 \"yellow\")\n");
	    	break;
	case -1:string com;
		if (mode==0) com="tot";
		else if (mode==1) com="cor";
		else if (mode==2) com="rel";
		fprintf ( Fgplot,	"unset xtics\n"
					"unset ytics\n"
					"unset border\n"
					"set colorbox user\n"
					"set style line 3432 lc \"white\"\n"
					"set colorbox border 3432\n"
					"set border ls 3432\n"
					"set colorbox vertical\n"
					"set colorbox size 0.05,0.20\n"
					"set colorbox origin 0.36,0.50\n"
					"set label \"%s. counts per \\n%dx%d channels\"  textcolor \"white\" at screen 0.285,screen 0.47 rotate front\n"
					"set label \"alpha channel\"  textcolor \"white\" at screen 0.6,screen 0.74-0.02 front\n"
					"set label \"gamma channel\"  textcolor \"white\" at screen 0.26+0.02,screen 0.1 rotate front\n"
					"replot\n",com.c_str(),step_alpha,step_gamma);
	    				break;
    	}
	fflush(Fgplot);
}

double minnum=100;
template <class inttype>
double _integrate(inttype *time_array, int mode, double dbkgnd, int subr)
{
	long long int result=0;
	long long int bkgnd=0;

	if (mode==1 || mode==2){
		for (unsigned i=0;i!=time_len/4;i++)          bkgnd+=time_array[i];
		for (unsigned i=3*time_len/4;i!=time_len;i++) bkgnd+=time_array[i];
		bkgnd*=2;
	}
	
	for (unsigned i=0;i!=time_len;i++)   
	        result+=time_array[i];
	        
	if (mode>=3){
		result=0;
		for (unsigned i=time_len/4;i!=3*time_len/4;i++)   
			result+=time_array[i];
	}
	
	if (subr==1) return result;
	if (mode==2){
		if (subr==2) return bkgnd;
		else if ((result-bkgnd)>=minnum) return ((double)(result-bkgnd)/bkgnd);
		else return NAN;
	}
	else if (mode==4){
		if (subr==2) return dbkgnd;
		else if ((result-dbkgnd)>=minnum) return ((double)(result-dbkgnd)/dbkgnd);
		else return NAN;
	}
	else if (mode==0 || mode==1) return (double)(result-bkgnd);
	else if (mode==3) return (double)(result-dbkgnd);
}

int main(int argc,char *argv[])
{
	if (argc!=3){_help();return 0;}
	int mode = atoi(argv[2]);
	if ((mode<0)||(mode>4)){_help();return 0;}
	string fpath = argv[1];
	if (fpath.back()=='/') fpath.resize(fpath.size()-1);
	FILE *Falpha,*Fgamma,*Ftime,*Fconf;
	Falpha = fopen((fpath + "/alpha.dat").c_str(),"rb");
	if (Falpha==NULL) {printf("Missing alpha.dat. Aborting...\n");exit(0);}
	Fgamma = fopen((fpath + "/gamma.dat").c_str(),"rb");
	if (Fgamma==NULL) {printf("Missing gamma.dat. Aborting...\n");exit(0);}
	Ftime = fopen((fpath + "/time.dat").c_str(),"rb");
	if (Ftime==NULL) {printf("Missing time.dat. Aborting...\n");exit(0);}
	_load_conf(fpath);
	
	unsigned *alpha_array = new unsigned[alpha_Nmes];
	unsigned *gamma_array = new unsigned[gamma_Nmes];
	fread (alpha_array,sizeof(unsigned),alpha_Nmes,Falpha);
	for (unsigned i=0;i!=alpha_Nmes-1*step_alpha;i++) if (alpha_max<alpha_array[i]) alpha_max=alpha_array[i];
	fclose(Falpha);
	fread (gamma_array,sizeof(unsigned),gamma_Nmes,Fgamma);
	for (unsigned i=0;i!=gamma_Nmes-1*step_gamma;i++) if (gamma_max<gamma_array[i]) gamma_max=gamma_array[i];
	fclose(Fgamma);
	
	long long unsigned alpha_tot=0;
	unsigned *alpha_array_sh = new unsigned[alpha_Nch];
	for (unsigned i=0;i!=alpha_Nch;i++) alpha_array_sh[i]=0;
	for (unsigned i=0;i!=alpha_Nmes;i++) {
		alpha_array_sh[(int)(i/step_alpha)]+=alpha_array[i];
		alpha_tot+=alpha_array[i];
	}
	long long unsigned gamma_tot=0;
	unsigned *gamma_array_sh = new unsigned[gamma_Nch];
	for (unsigned i=0;i!=gamma_Nch;i++) gamma_array_sh[i]=0;
	for (unsigned i=0;i!=gamma_Nmes;i++) {
		gamma_array_sh[(int)(i/step_gamma)]+=gamma_array[i];
		gamma_tot+=gamma_array[i];
	}
	long long unsigned bkgnd_tot=0;
	unsigned *time_array = new unsigned[time_len];
	for (unsigned i=0;i!=alpha_Nch;i++){
		for (unsigned j=0;j!=gamma_Nch;j++){
			fread (time_array,sizeof(unsigned),time_len,Ftime);
			for (unsigned k=0;k!=time_len/4;k++)          bkgnd_tot+=time_array[k];
			for (unsigned k=3*time_len/4;k!=time_len;k++) bkgnd_tot+=time_array[k]; 
		}
	}
	fseek (Ftime, 0,SEEK_SET);
	
	bool done=false;
	bool done2=false;
	_gnuplot_set(0,mode);
plotgt:	_gnuplot_set(1,mode);
	fprintf ( Fgplot, "plot \"-\" using 1:2 with filledcurves below x lc rgb \"#19afe0\" title \"alpha\", \"-\" using 1:2 w l lc rgb \"#19afe0\" notitle\n");
	for (int k=0;k!=2;k++){	//do it 2x
		for (unsigned i=0;i!=alpha_Nmes;i++){
			if (i<=(alpha_Nch-1)*step_alpha) fprintf(Fgplot,"%u %u\n",i,alpha_array[i]);
		} 
		fprintf ( Fgplot, "e\n" );
	}
	fflush(Fgplot);
	
	_gnuplot_set(2,mode);
	fprintf ( Fgplot, "plot \"-\" using 1:2 with filledcurves below y2 lc rgb \"#74c22b\" notitle, \"-\" using 1:2 w l lc rgb \"#74c22b\" notitle\n");
	for (int k=0;k!=2;k++){	//do it 2x
		for (unsigned i=0;i!=gamma_Nmes;i++){
			if (i<=(gamma_Nch-1)*step_gamma) fprintf(Fgplot,"%u %u\n",gamma_array[i],i);
		} 
		fprintf ( Fgplot, "e\n" );
	}
	fflush(Fgplot);
	
	_gnuplot_set(-2,mode);
	fprintf ( Fgplot, "splot \"-\" using 1:2:3 with pm3d notitle\n");
	double sum;
	for (unsigned i=0;i!=alpha_Nch;i++){
		for (unsigned j=0;j!=gamma_Nch;j++){
			fread (time_array,sizeof(unsigned),time_len,Ftime);
			sum = _integrate<unsigned>(time_array, mode, (double)bkgnd_tot/alpha_tot*alpha_array_sh[i]/gamma_tot*gamma_array_sh[j],0);
			if (sum<0) sum=0;
			if (j<=gamma_Nch-1 && i<=alpha_Nch-1) fprintf(Fgplot,"%lf %lf %lf\n",(double)i+0.5,(double)j+0.5,sum);
		}
		if (i<=alpha_Nch-1)fprintf ( Fgplot, "\n" );
	} 
	fprintf ( Fgplot, "e\n" );
	fflush(Fgplot);
	_gnuplot_set(-1,mode);
	
	if (!done){
		char newvalstr[100]; double val;
		printf ( "redefine alpha min? (hover over colormap plot to get value)(input a value or press enter to keep current setting)\n");char * gets ( char * str );
		fgets (newvalstr , 100 , stdin);
		if (sscanf(newvalstr,"%lf",&val)==1) xmin=val-0.5;
		printf ( "redefine alpha max?\n");
		fgets (newvalstr , 100 , stdin);
		if (sscanf(newvalstr,"%lf",&val)==1) xmax=val-0.5;
		printf ( "redefine gamma min?\n");
		fgets (newvalstr , 100 , stdin);
		if (sscanf(newvalstr,"%lf",&val)==1) ymin=val-0.5;
		printf ( "redefine gamma max?\n");
		fgets (newvalstr , 100 , stdin);
		if (sscanf(newvalstr,"%lf",&val)==1) ymax=val-0.5;
		printf ( "redefine alpha tics step? (IN CHANNELS - look at the other two graphs)\n");
		fgets (newvalstr , 100 , stdin);
		if (sscanf(newvalstr,"%lf",&val)==1) xstep=val;
		printf ( "redefine gamma tics step?\n");
		fgets (newvalstr , 100 , stdin);
		if (sscanf(newvalstr,"%lf",&val)==1) ystep=val;
		double cbtics=-1;
		printf ( "redefine colorbox tics step?\n");
		fgets (newvalstr , 100 , stdin);
		if (sscanf(newvalstr,"%lf",&val)==1) cbtics=val;
		
		done=true;
		fprintf ( Fgplot,"reset\nclear\n");
		if (cbtics!=-1) fprintf ( Fgplot,"set cbtics %lf\n",cbtics);
		fflush(Fgplot);
		
		fseek (Ftime, 0,SEEK_SET);
		goto plotgt;
	}
	
	
	
	if (mode==2 || mode==4){
		if (!done2) {
			printf ( "\n\nEnter threshold and tics step and hit enter, repeat until satisfied . To end hit enter again.\n");
			done2=true;
		}
		char newvalstr[100]; double cbtics;
		fgets (newvalstr , 100 , stdin);
		if (sscanf(newvalstr,"%lf%lf",&minnum,&cbtics)==2){
			fseek (Ftime, 0,SEEK_SET);
			fprintf ( Fgplot,"reset\nclear\n");
			fprintf ( Fgplot,"set cbtics %lf\n",cbtics);
			fflush(Fgplot);
			goto plotgt;
		}
	}
	else {
		unsigned *exptime_array = new unsigned[time_len];
		printf ( "\n\nSet areas by providing boundaries. Use format: amin amax gmin gmax and hit enter. To end hit enter again.\n");
		for (int i=1;;i++){
			char newvalstr[100]; double amin,amax,gmin,gmax;
			printf ( "Area%d:\n",i);
			fgets (newvalstr , 100 , stdin);
			if (sscanf(newvalstr,"%lf%lf%lf%lf",&amin,&amax,&gmin,&gmax)!=4) break;
			fprintf ( Fgplot,	"set label \"Area%d\" at %lf,%lf textcolor \"white\" front\n"
						"set arrow from %lf,%lf to %lf,%lf nohead dt 3 lc \"white\" front\n"
						"set arrow from %lf,%lf to %lf,%lf nohead dt 3 lc \"white\" front\n"
						"set arrow from %lf,%lf to %lf,%lf nohead dt 3 lc \"white\" front\n"
						"set arrow from %lf,%lf to %lf,%lf nohead dt 3 lc \"white\" front\n"
						"replot\n"
						,i,amax+(xmax-xmin)/100,gmax-(ymax-ymin)/40,
						amin,gmin,amax,gmin,
						amax,gmin,amax,gmax,
						amax,gmax,amin,gmax,
						amin,gmax,amin,gmin);
			fflush(Fgplot);	
			
			double *expalpha_array = new double[(int)(amax-amin+2)];
			for (int j=0;j<amax-amin+1;j++) expalpha_array[j]=0;
			double *expgamma_array = new double[(int)(gmax-gmin+2)];
			for (int j=0;j<gmax-gmin+1;j++) expgamma_array[j]=0;
			
			long unsigned *expalpha_array_time;
			expalpha_array_time = new long unsigned[time_len];
			
			long unsigned **expgamma_array_time;
			expgamma_array_time = new long unsigned *[(int)(gmax-gmin+2)];
			for (int j=0;j<gmax-gmin+1;j++){
				expgamma_array_time[j] = new long unsigned[time_len];
				for (int k=0;k!=time_len;k++) expgamma_array_time[j][k]=0;
			} 
			
			long unsigned sum;
			
			fseek (Ftime, 0,SEEK_SET);
			long long unsigned total=0;
			for (unsigned i=0;i!=time_len;i++) exptime_array[i]=0;

			for (unsigned i=0;i!=alpha_Nch;i++){
				sum=0;
				for (int j=0;j!=time_len;j++) expalpha_array_time[j]=0;
				for (unsigned j=0;j!=gamma_Nch;j++){
					fread (time_array,sizeof(unsigned),time_len,Ftime);
					if ((i>=amin)&&(i<=amax)&&(j>=gmin)&&(j<=gmax)){
						for (unsigned k=0;k!=time_len;k++) {
							exptime_array[k]+=time_array[k];
							total+=time_array[k];
							expalpha_array_time[k]+=time_array[k];
							expgamma_array_time[j-(int)gmin][k]+=time_array[k];
						}
					}
					sum+=gamma_array_sh[j];
				}
				if ((i>=amin)&&(i<=amax)) expalpha_array[i-(int)amin]=_integrate<long unsigned>(expalpha_array_time, mode, (double)bkgnd_tot/alpha_tot*alpha_array_sh[i]/gamma_tot*sum,0); 
			}

			for (unsigned j=0;j!=gamma_Nch;j++){
				sum=0;
				for (unsigned i=0;i!=alpha_Nch;i++) sum+=alpha_array_sh[i];
				if ((j>=gmin)&&(j<=gmax)) expgamma_array[j-(int)gmin]=_integrate<long unsigned>(expgamma_array_time[j-(int)gmin], mode, (double)bkgnd_tot/alpha_tot*sum/gamma_tot*gamma_array_sh[j],0);
			}
			FILE *Ftimesum;
			string fname="time_area"+to_string(i)+".dat";
			Ftimesum = fopen(fname.c_str(),"wb");
			fwrite (exptime_array,sizeof(unsigned),time_len,Ftimesum);
			fclose(Ftimesum);
			printf ( "Total number of counts in area: %llu\n",total);
			printf ( "Time graph exported to file %s\n",fname.c_str());

			fname="alpha_area"+to_string(i)+".txt";
			Ftimesum = fopen(fname.c_str(),"w");
			for (int j=0;j<amax-amin+1;j++) fprintf(Ftimesum,"%lf %lf\n",(amin+j+0.5)*step_alpha,expalpha_array[j]);
			fclose(Ftimesum);
			printf ( "Alpha graph exported to file %s\n Format is: Channel CountN\n",fname.c_str());
			delete[] expalpha_array;
			delete[] expalpha_array_time;

			fname="gamma_area"+to_string(i)+".txt";
			Ftimesum = fopen(fname.c_str(),"w");
			for (int j=0;j<gmax-gmin+1;j++) fprintf(Ftimesum,"%lf %lf\n",(gmin+j+0.5)*step_gamma,expgamma_array[j]);
			fclose(Ftimesum);
			printf ( "Gamma graph exported to file %s\n Format is: Channel CountN\n",fname.c_str());
			delete[] expgamma_array;
			for (int j=0;j<gmax-gmin+1;j++) delete[] expgamma_array_time[j];
			delete[] expgamma_array_time;
		}
	}
		
	
	
	fclose(Ftime);

	printf ( "\n\npress enter to exit (gnuplot window will lose function)\n");
	scanf ("%*c");
	fclose(Fgplot);
}
