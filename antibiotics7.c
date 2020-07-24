#include <stdio.h>
#include <string.h>
#include <time.h>
#include <errno.h>
#include <math.h>
#include <grace_np.h>
#include <unistd.h>
#include <float.h>
#include <limits.h>
#include <signal.h>
#include <cash2003.h>
#include <cash2.h>
#include <mersenne.h>
#include <cash2-s.h>

static TYPE2** Medium;
static TYPE2** Cost;
static TYPE2** ColorMap;
static TYPE2** Diffusion_plane;

static TYPE2 empty={0,0,0,0,0,0.,0.,0.,0.,0.};



double growthrate=0.6; // Speed of growth, chance to invade empty space
float death=0.07; // Death rate
float minimum_cost=0.1;
double initial_cost = 0.1; // Cost of antibiotic production
float mut_step_cost=0.1;
float mut_rate_cost=0.3;
float mut_step_production=0.1;
float mut_rate_production=0.3;
float diffusion=0.2;
double decayrate = 0.03;
double Rearrangement_rate=0.015;
double growthrate_competitor = 0.6;
double initial_production = 0.3;
int killing_time = 10000;
int number_sensitives = 1000;
int killing_period = 0;
double nastiness = 2.0;
int wt_seeds = 15;



char *name;
char buf[400]; //<-- temp buffer for string
char folder[400]; //<-- another temp buffer for string
char gridfolder[400]; //<-- another temp buffer for string
char costsfolder[400];
char productionfolder[400];
// char name[30] = "data";
void Initial(void)
{
  display = 0;
  char* readOut; // Stores arguments for loop
  for(int i = 0; i < (int)argc_g; i++) // Loops over the words on the command line
  {
    readOut = (char*)argv_g[i]; // Stores current word
    if(strcmp(readOut, "-mc") == 0) minimum_cost = atof(argv_g[i+1]);            // atof voor charater-to-float (of double)
    if(strcmp(readOut, "-n") == 0) nastiness = atof(argv_g[i+1]);         // atoi voor character-to-integer
    if(strcmp(readOut, "-d") == 0) diffusion = atof(argv_g[i+1]);
    if(strcmp(readOut, "-dec") == 0) decayrate = atof(argv_g[i+1]);
    if(strcmp(readOut, "-X") == 0) display = 0;
    if(strcmp(readOut, "-rr") == 0) Rearrangement_rate = atof(argv_g[i+1]);
    if(strcmp(readOut, "-ic") == 0) initial_cost = atof(argv_g[i+1]);
    if(strcmp(readOut, "-ip") == 0) initial_production = atof(argv_g[i+1]);
    if(strcmp(readOut, "-name") == 0) name = argv_g[i+1];
    if(strcmp(readOut, "-r") == 0) ulseedG = atoi(argv_g[i+1]);
    if(strcmp(readOut, "-gr") == 0) growthrate = atof(argv_g[i+1]);
    if(strcmp(readOut, "-gc") == 0) growthrate_competitor = atof(argv_g[i+1]);
    if(strcmp(readOut, "-ns") == 0) number_sensitives = atoi(argv_g[i+1]);
  }

  sprintf(folder,"/linuxhome/tmp/bramve/%s/",name);
  sprintf(buf,"mkdir %s -p",folder);
  int success = system(buf);
  if(success == -1) {
    printf("Error in system call to make dir...\n"); exit(0);
  } //Exits if fails.

  init_genrand(ulseedG); // random seed moet weer geinitialiseerd worden volgens Bram

  sprintf(gridfolder,"/linuxhome/tmp/bramve/%s/grid",name);
  sprintf(buf,"mkdir %s -p",gridfolder);
  success = system(buf);
  if(success == -1) {
    printf("Error in system call to make grid dir...\n"); exit(0);
  } //Exits if fails.

  sprintf(productionfolder,"/linuxhome/tmp/bramve/%s/grid",name);
  sprintf(buf,"mkdir %s -p",gridfolder);
  success = system(buf);
  if(success == -1) {
    printf("Error in system call to make grid dir...\n"); exit(0);
  } //Exits if fails.

  sprintf(costsfolder,"/linuxhome/tmp/bramve/%s/grid",name);
  sprintf(buf,"mkdir %s -p", costsfolder);
  success = system(buf);
  if(success == -1){
    printf("Error in system call to make costs dir...\n"); exit(0);}

	MaxTime = 96000000;
	nrow = 400;
	ncol = 400;
	nplane =4;
	scale = 1;
	boundary = WRAP;
	//ulseedG = 57;

	boundaryvalue2 = (TYPE2){0,0,0,0,0,0.,0.,0.,0.,0.};
}

void Restart_Simulation(char* filename)
{
  FILE* datafile;
  char* Mediumval;
  char* Mediumfval;
  char* Mediumfval2;
  char* Diffusionfval;
  char* x;
  char* y;
  char line[1000];
  int i, j;
  datafile = fopen(filename, "r");

  if(datafile == NULL){
    fprintf(stderr, "Error opening file\n");
    exit(1);
  }
  for(i=1;i<=nrow;i++){
    for(j=1;j<=ncol;j++){
      fgets(line, 1000, datafile);
      Mediumval = strtok(line,"\t");
      Mediumfval = strtok(NULL, "\t");
      Mediumfval2 = strtok(NULL, "\t");
      Diffusionfval = strtok(NULL, "\t");
      Mediumfval3 = strtok(NULL, "\t")
      x = strtok(NULL, "\t");
      y = strtok(NULL, "\t");
      int row = atoi(x);
      int col = atoi(y);
      Medium[row][col].val = atoi(Mediumval);
      Medium[row][col].fval = atof(Mediumfval);
      Medium[row][col].fval2 = atof(Mediumfval2);
      Diffusion_plane[row][col].fval = atof(Diffusionfval);
    }
  }
  fclose(datafile);
}

void InitialPlane(void)
{


	MakePlane(&Medium, &Cost, &Diffusion_plane, &ColorMap);
		TYPE2 bg, s1, s2, s3;
		bg.val = 0;
		s1.val = 1;
		s2.val = 2;
    s3.val = 3;
		bg.fval = 0.;
		s1.fval = initial_cost;
    s1.fval2 = initial_production;
		s2.fval = 0.;
    s2.fval2 = 0.;
    s3.fval = 0.;
    s3.fval2 = 0.;
    double chance = 0.0000625; // very scientific
    printf("%f", chance);
		InitialSetS(Medium, 0, bg);
    InitialSet(Diffusion_plane, 0, 0);
		ColorRGB(3,0,0,255);

	int m,n;
	for(m=1;m<=nrow;m++)
	for(n=1;n<=ncol;n++){
		ColorMap[m][n].val = 110-(float)m/(float)nrow*100;
	}

	int r=0,g=0,b=255;
  double nr=102.;    //nr of ColorRGB's you want to create
  double range=1275.;  // range of coloursteps: 255*5= 1275
  double x = range/nr;
  int y=0,c;
	for(c=0;c<(int)range;c++){
	 if(c<255){			//starts blue
		 r = r + 1;		//goes magenta
	 }else if(c<255*2){
		 b = b - 1;		//goes red
	 }else if(c<255*3){
		 g = g + 1;		//goes yellow
	 }else if(c<255*4){
		 r = r -1;    		//goes green
	 }else if(c<255*5){
		 b = b + 1;		//goes cyan
	 }

	 if(c == (int)(y*x+0.5)){
		 ColorRGB(10+y,r,g,b);	//so colour indexes in the gradient start at 10
		 y++;
	 }
 }

}


 void NextState(int row, int col)
 {
   TYPE2 neighbor = empty;
   neighbor = RandomMooreS8(Medium, row, col);

   if(Medium[row][col].val == 0 && neighbor.val == 1){
     if(genrand_real1() < growthrate){
       if(genrand_real1() > Rearrangement_rate){ // WT reproduction
         Medium[row][col] = neighbor;
         if(genrand_real1() < mut_rate_cost){ // WT mutates cost parameter
           Medium[row][col].fval += mut_step_cost * (genrand_real2()/10 - 0.05);
           if(Medium[row][col].fval < minimum_cost){ // bounding cost parameter
             Medium[row][col].fval = minimum_cost;
           }
         }
         if(genrand_real1() < mut_rate_production){ // WT mutates production parameter
           Medium[row][col].fval2 += mut_step_production * (genrand_real2()/10 - 0.05);
         }
       }


			else {
				Medium[row][col].val = 2 ;
				Medium[row][col].fval = neighbor.fval;
        Medium[row][col].fval2 = neighbor.fval2;
        if(genrand_real1() < mut_rate_cost){ // cost parameter mutates
          Medium[row][col].fval += mut_step_cost * (genrand_real2()/10 - 0.05);
            if(Medium[row][col].fval < minimum_cost){
              Medium[row][col].fval = minimum_cost;
            }
        }
        if(genrand_real1() < mut_rate_production){ // mutant mutates production parameter
          Medium[row][col].fval2 += mut_step_production * (genrand_real2()/10 - 0.05);
        }
      }
    }
  }
  // Mutant reproduction
  if(Medium[row][col].val == 0 && neighbor.val == 2){
  	if(genrand_real1() < (growthrate - neighbor.fval)){
  		Medium[row][col] = neighbor;
      // mutations
      if(genrand_real1() < mut_rate_cost){
        Medium[row][col].fval += mut_step_cost * (genrand_real2()/10 - 0.05);
        if(Medium[row][col].fval < minimum_cost)
          Medium[row][col].fval = minimum_cost;
      }
      if(genrand_real1() < mut_rate_production){ // mutant mutates nastiness parameter
        Medium[row][col].fval2 += mut_step_production* (genrand_real2()/10 - 0.05);
      }
    }
  }


  if(Medium[row][col].val == 1 || Medium[row][col].val == 2){

  	if(genrand_real1() < death){
  		Medium[row][col] = empty;
    }


    if(Medium[row][col].val == 2){

  		Diffusion_plane[row][col].fval += Medium[row][col].fval2;
    }
  }


if(Medium[row][col].val == 0 && neighbor.val == 3){
  if(genrand_real1() < growthrate){
    Medium[row][col].val = 3;}}

if(Medium[row][col].val == 3){
  if(genrand_real1() < death + nastiness*Diffusion_plane[row][col].fval){
    Medium[row][col].val = 0;}
  }

if(Diffusion_plane[row][col].fval > 0.0){
	Diffusion_plane[row][col].fval = Diffusion_plane[row][col].fval *(1 - decayrate);
}

}

int GetColorIndexFrom(int val, double fval, double scale_max)
{
		int color;
		double max;
		max = scale_max;
		if(fval > max){
			color = 100;}
		else color = (int)(100.*fval/max)+10;
		return color;
}

void Update(void)
{
  int i, j;
  DiffusionFVAL(Diffusion_plane,diffusion,1);

  //Killing
  if(Time >0 && Time%killing_time == 0){
    PerfectMix(Medium);


    int counter = 0;
    int arr_x[15];
    int arr_y[15];
    int arr_type[15];
    double arr_cost[15];
    double arr_production[15];
    while(counter < wt_seeds){
      int x = genrand_int(1, nrow);
      int y = genrand_int(1, nrow);
      if(Medium[x][y].val == 1 || Medium[x][y].val == 2){
        arr_x[counter] = x;
        arr_y[counter] = y;
        arr_cost[counter] = Medium[x][y].fval;
        arr_production[counter] = Medium[x][y].fval2;
        arr_type[counter] = Medium[x][y].val;
        counter +=1;
      }
    }
    int i, j;
		for(i=1;i<=nrow;i++)
		for(j=1;j<=ncol;j++){
      Diffusion_plane[i][j].fval = 0;
      if(Medium[i][j].val == 3) Medium[i][j] = empty;
      if(Medium[i][j].val == 2) Medium[i][j] = empty;
      if(Medium[i][j].val == 1) Medium[i][j] = empty;
    }
    for(i=0;i<=wt_seeds-1;i++){
      if(arr_type[i] == 1){
        Medium[arr_x[i]][arr_y[i]].val = 1;
        Medium[arr_x[i]][arr_y[i]].fval = arr_cost[i];
        Medium[arr_x[i]][arr_y[i]].fval2 = arr_production[i];
      }
    }




      //if(genrand_real1() > 0.000425) Medium[i][j] = empty;

    PerfectMix(Medium);
    int additor = 0;
    while(additor < number_sensitives){
      int i = genrand_int(1, nrow);
      int j = genrand_int(1, ncol);
      if(Medium[i][j].val == 0){
        Medium[i][j].val = 3;
        Medium[i][j].fval = 0;
        additor++;
         }
	     }
     }


  FILE *restartfile;
  if(Time%1000000){
    char restartbuffer[400];
    sprintf(restartbuffer,"%s/restart.txt",gridfolder);
    restartfile = fopen(restartbuffer, "a");
    for(i=1;i<=nrow;i++)
    {
      for(j=1;j<=ncol;j++)
      {
        fprintf(restartfile, "%d\t%f\t%f\t%f\t%d\t%d\n", Medium[i][j].val, Medium[i][j].fval, Medium[i][j].fval2, Diffusion_plane[i][j].fval, i, j);
      }
    }
    fclose(restartfile);
  }

  FILE *populationfile;
  if(Time%2500 == 0){
    char populationbuffer[400];
    sprintf(populationbuffer,"%s/popsize2.txt",folder);    //T=1000 "/linuxhome/tmp/brem/Hierzo/grid_at_time_1000.dat"
    populationfile = fopen(populationbuffer, "a");
    if(Time == 0){
      fprintf(populationfile, "Time\tPopulation_wildtype\tPopulation_mutants\tsensitives");
      fprintf(populationfile, "\n");}
    int pop_WT = countGlobal(Medium, 1);
    int pop_mutant = countGlobal(Medium, 2);
    int pop_sensitive = countGlobal(Medium, 3);
    fprintf(populationfile, "%d\t%d\t%d\t%d\n", Time, pop_WT, pop_mutant, pop_sensitive);
    fclose(populationfile);
  }

  FILE *costsfile;

  if(Time%1000 == 0){
    char costsbuffer[400];
    sprintf(costsbuffer,"%s/costs_concentration_production.txt",folder);    //T=1000 "/linuxhome/tmp/brem/Hierzo/grid_at_time_1000.dat"
    costsfile = fopen(costsbuffer, "a");
    if(Time == 0){
      fprintf(costsfile, "Time\tCost\tConcentration\tProduction\tI\tJ\n");
      }
    for(int r=1;r<=120;){
       i = genrand_int(1,nrow);
       j = genrand_int(1,ncol);
       if(Medium[i][j].val > 0 && Medium[i][j].val != 3){
         fprintf(costsfile, "%d\t%f\t%f\t%f\t%d\t%d\n", Time, Medium[i][j].fval,Diffusion_plane[i][j].fval, Medium[i][j].fval2,i,j);
         r +=1;
       }
     }
    fclose(costsfile);
  }

  if(Time%killing_time == 0){
    killing_period += killing_time;
  }

    FILE *valfile;
    if(Time%2500 == 0)
    {
      char valbuffer[400];
      sprintf(valbuffer,"%s/grid_type%d.txt", gridfolder, Time);
      valfile = fopen(valbuffer, "a");
      for(i=1;i<=nrow;i++)
      {
        for(j=1;j<=ncol;j++)
        {
          fprintf(valfile, "%d\t%d\t%d\n", Medium[i][j].val, i, j);
        }
      }
      fprintf(valfile, "\n");
      fclose(valfile);
    }

      FILE *evolutionfile;
      if(Time%2500 == 0)
      {
        char evolutionbuffer[400];
        sprintf(evolutionbuffer,"%s/grid_costs%d.txt", costsfolder, Time);
        evolutionfile = fopen(evolutionbuffer, "a");
        for(i=1;i<=nrow;i++)
        {
          for(j=1;j<=ncol;j++)
          {
            fprintf(evolutionfile, "%f\t%d\t%d\n", Medium[i][j].fval, i, j);
          }
        }
        fprintf(evolutionfile, "\n");
        fclose(evolutionfile);
      }

      FILE *productionfile;
      if(Time%2500 == 0)
      {
        char productionbuffer[400];
        sprintf(productionbuffer,"%s/grid_production%d.txt", costsfolder, Time);
        productionfile = fopen(productionbuffer, "a");
        for(i=1;i<=nrow;i++)
        {
          for(j=1;j<=ncol;j++)
          {
            fprintf(productionfile, "%f\t%d\t%d\n", Medium[i][j].fval2, i, j);
          }
        }
        fprintf(productionfile, "\n");
        fclose(productionfile);
      }

    FILE *runfile;
    char commandfile[400];
    sprintf(commandfile,"%s/command%s.txt",folder, name);
    runfile = fopen(commandfile, "w");
    for(int i = 0; i < (int)argc_g; i++){
      fprintf(runfile,"%s ", (char*)argv_g[i]);}
    fclose(runfile);

    if(Time%50 ==0)
    {

      Synchronous(1, Medium);
      int i,j;
      for(i=1;i<=nrow;i++) for(j=1;j<=ncol;j++)
      {
        Diffusion_plane[i][j].val=GetColorIndexFrom(Diffusion_plane[i][j].val, Diffusion_plane[i][j].fval, 2.);
        Cost[i][j].val=GetColorIndexFrom(Medium[i][j].val, Medium[i][j].fval, 0.6);
      }
      Display(Medium, Diffusion_plane, Cost, ColorMap);
    }

}
