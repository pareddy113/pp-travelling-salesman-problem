/*
 The parallel Branch and bound is implemented as follows
  Nodes are assigned such that each processor gets one node each
  If size/nprocs is > 1, then each processor is assigned size/nprocs number of nodes
  Each node is taken as the current starting node
  DFS is done on this node to find its nearest neighbor based on the cost metric.
 Further searching of the current branch are pruned by using the cost of the best
 solution so far.
  If the cost of traveling to this city is less than the global best, the new node is taken
 as the current node and DFS is performed on this node
  This process goes on till all the nodes have been visited and we visit the root node
 back
 
 */

#include <float.h>
#include <getopt.h>
#include <limits.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

int n, lrank, nprocs;
int compstat = 0;
char locsf[255];
double gen_time, proc_time, comm_time, total_time;

typedef enum
 {
   linear,
   nn
 } TYPE;
TYPE type;

typedef struct
 {
   int line, v;
   float x, y;
 } LOCATION;
LOCATION *locs;

typedef struct
 {
   int pointA, pointB;
   float cost;
 } COST;

void swap(int *p1, int *p2);
float calcost(int *a);
void display(int *a, int *c, float cost);
float distance(LOCATION a, LOCATION b);
void permute();
int nearest(int curr, int start, int end);
void n_neigh();
int par_args(int argc, char **argv);
void par_file();

int main(int argc, char **argv)
 {
   //printf("start\n ");
   int i;
   double t_start, t_end;
   float gen_time = 0.0, proc_time = 0.0, comm_time = 0.0, total_time = 0.0;
   if(par_args(argc, argv))
      return 1;
  // printf("start process\n");
  // return 0;
   par_file();
   //return 0;
  // printf("read file\n");
   MPI_Init(&argc, &argv);
   MPI_Comm_rank(MPI_COMM_WORLD, &lrank);
   MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
   t_start = MPI_Wtime();
   if(type == linear)
      permute();
   else if(type == nn)
      n_neigh();
   t_end = MPI_Wtime();
   total_time = t_end - t_start;
  // printf("finished processing\n");
   if(compstat)
    {
      printf("%d\tg\t%d\t%d\t%f\n", n, lrank, nprocs, gen_time);
      printf("%d\tp\t%d\t%d\t%f\n", n, lrank, nprocs, proc_time);
      printf("%d\tc\t%d\t%d\t%f\n", n, lrank, nprocs, comm_time);
      printf("%d\tt\t%d\t%d\t%f\n", n, lrank, nprocs, total_time);
    }
   free(locs);
   MPI_Finalize();
   return 0;
 }

void swap(int *p1, int *p2)
 {
   int temp;
   temp = *p1;
   *p1 = *p2;
   *p2 = temp;
 }
float calcost(int *a)
 {
   int i;
   float cost = 0.0f;
   for(i = 0; i < n - 1; i++)
       cost += distance(locs[a[i] - 1], locs[a[i + 1] - 1]);
   return cost;
 }

void display(int *a, int *c, float cost)
 {
   int x;
   for(x = 0; x < n; x++)
       printf("%d  ",a[x]);
   printf("cost:%f c:%d\n", cost, *c);
 }

float distance(LOCATION a, LOCATION b)
 {
   return sqrt(pow(a.x - b.x, 2) + pow(a.y - b.y, 2));
 }

void permute()
 {
   int i, c, x, y , np = 1, *a;
   double start, end, dt;
   float cost = 0.0;
   start = MPI_Wtime();
   for(i = 1; i <= n; i++)
       np *= i;
   a = (int*) malloc(sizeof(int) * n * np);
   for(i = 0; i < n; i++)
       a[i] = i + 1;
   while(c < np)
    {
       for(y = 0; y < n - 1; y++)
        {
           swap(&a[y], &a[y + 1]);
           cost = calcost(a);
           if(!compstat)
              display(a, &c, cost);
           c++;
        }
       swap(&a[0], &a[1]);
       cost = calcost(a);
       if(!compstat)
          display(a, &c, cost);
       c++;
       for(y = n - 1; y > 0; y--)
        {
           swap(&a[y], &a[y - 1]);
           cost = calcost(a);
           if(!compstat)
              display(a, &c, cost);
           c++;
        }
       swap(&a[n - 1], &a[n - 2]);
       cost = calcost(a);
       if(!compstat)
          display(a, &c, cost);
          c++;
    }
   end = MPI_Wtime();
   dt = end - start;
   proc_time += dt;
 }

int nearest(int curr, int start, int end)
 {
   int i, index = -1;
   float min = FLT_MAX;
   float dist;
   for(i = start; i<= end; ++i)
    {
       dist = distance(locs[curr], locs[i]);
       if(dist < min && i != curr && locs[i].v == 0)
        {
          min = dist;
          index = i;
        }
    }
   return index;
 }
void n_neigh()
 {
   int i, j, index, sloc, eloc, next;
   int locpn = n / nprocs;
   int *inm;
   int fpath[n];
   double start, end, dt;
   float min = FLT_MAX;
   float dist;
   float cost = 0.0f;
   inm = (int*)malloc(sizeof(int) * nprocs);
   sloc = locpn * lrank;
   eloc = sloc + locpn - 1;
   if(lrank == nprocs - 1)
      eloc += n % nprocs;
   next = 0;
   fpath[0] = 0;
   for(i = 0; i < n - 1; i++)
    {
       start = MPI_Wtime();
       MPI_Bcast(&next, 1, MPI_INT, 0, MPI_COMM_WORLD);
       end = MPI_Wtime();
       dt = end - start;
       comm_time += dt;
       start = MPI_Wtime();
       locs[next].v = 1;
       int index = nearest(next, sloc, eloc);
       end = MPI_Wtime();
       dt = end - start;
       proc_time += dt;
       start = MPI_Wtime();
       MPI_Gather(&index, 1, MPI_INT, inm, 1, MPI_INT, 0, MPI_COMM_WORLD);
       end = MPI_Wtime();
       dt = end - start;
       comm_time += dt;
       if(lrank == 0)
        {
          start = MPI_Wtime();
          index = inm[0];
          min = FLT_MAX;
          for(j = 0; j < nprocs; ++j)
           {
              if(inm[j] < 0)
                 continue;
              dist = distance(locs[next], locs[inm[j]]);
              if(dist < min)
               {
                 min = dist;
                 index = inm[j];
               }
           }
          next = index;
          fpath[i + 1] = index;
          end = MPI_Wtime();
          dt = end - start;
          proc_time += dt;
        }
       MPI_Barrier(MPI_COMM_WORLD);
    }
   if(lrank == 0 && !compstat)
    {
      for(i = 0; i < n; ++i)
          printf("%d ", fpath[i]);
      printf("\n");
    }
   free(inm);
 }

int par_args(int argc, char **argv)
 {
    int i, c, option_index = 0;
    char *result = NULL;
    char delims[] = "m";
    static struct option long_options[] =
     {
         {"input_method1", required_argument,  0, 'q'},
         {"input_method2", required_argument,  0, 'w'},
         {0, 0, 0, 0}
     };
    while((c = getopt_long (argc, argv, "f:n:t:c", long_options, &option_index)) != -1)
     {    
//	   printf("inside switch\n");
           switch(c)
            {
                  case 'c': compstat = 1;
//                            printf("case c compstat %d\n",compstat);
  			    break;
                  case 'f': strcpy(locsf, optarg);
  //                          printf("case f locsf %s\n",locsf);
                            break;
                            
                  case 'n': n = atoi(optarg);
    //                        printf("case n %d\n",n);
                            break;
                  case 't': if(strcmp(optarg, "linear") == 0)
                               type = linear;
                            else if(strcmp(optarg, "nn") == 0)
//{                            printf("case t %d \n",type);
                               type = nn;//}
                            else
                             {
                               fprintf(stderr, "Option -%c %s in incorrect. Allowed values are: linear, nn\n", optopt, optarg);
                               return 1;
                             }
//                           printf("case t %d\n",type);
                            break;
                  case '?': if(optopt == 'n')
                               fprintf(stderr, "Option -%c requires an argument.\n", optopt);
                            else if(isprint (optopt))
                               fprintf(stderr, "Unknown option `-%c'.\n", optopt);
                            else
                               fprintf (stderr, "Unknown option character `\\x%x'.\n", optopt);
                            return 1;
                  default:  fprintf(stderr, "Usage: %s -n <number of numbers> \n", argv[0]);
                            fprintf(stderr, "\tExample: %s -n 1000\n", argv[0]);
                            return 1;
            }
     }
printf("\n Processing Time %f \t Communication Time %f \n", proc_time, comm_time);
//printf("%d ,%d, %s, %d\n",compstat,n,locsf,type);
 return 0;
 }

void par_file()
 {
//printf("%d ,%d, %s, %d\n",compstat,n,locsf,type);
 
   FILE *fp;
   int i, line;
   char buff[1024];
   float x, y;
   fp = fopen(locsf, "r");
   for(i = 0; i < 9; i++)
       fgets(buff, 1024, fp);
  // printf("read file\n");
   while(fscanf(fp, "%d %f %f", &line, &x, &y) > 0 )
    {
         if(line == n)
                                   break;
    } 
//	printf("scanned file\n");
  locs = (LOCATION *) malloc(sizeof(LOCATION) * line);
   rewind(fp);
  // printf("rewind \n");
   for(i = 0; i < 9; i++)
       fgets(buff, 1024, fp);
	//printf("ignore 7 lines %d %s \n", line,buff);

    while(fscanf(fp, "%d %f %f", &line, &x, &y) > 0  && line <= n)
     {
          locs[line - 1].line = line;
          locs[line - 1].x = x;
          locs[line - 1].y = y;
          locs[line - 1].v = 0;
     }
    fclose(fp);

 }

