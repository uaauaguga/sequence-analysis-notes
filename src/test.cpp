#include<stdio.h>
#include<stdlib.h>
#include<getopt.h>
#include<string.h>
int main(int argc,char * argv[]){
int opt,age;
char *name;
/* 
 * val: the value to return is flag is NULL
 * the calling program may set val to the equivalent short option character, 
 * to make it compatible with short arguments in the while - switch case loop
 * Note the value of current parameter is always stored in the external variable optparg
*/
static option long_options[] = {
	{"name",required_argument,NULL,'n'},
	{"age", required_argument,NULL,'a'},
    {"help",no_argument,NULL,'h'}
};
int option_index = 0;
//printf("short-opt\tlong-opt-index\tcmd index\tvalue\tlong opt\tlong value\n");
while((opt=getopt_long(argc,argv,"n:a:h",long_options,&option_index))!=-1){
	//printf("%c\t%d\t%d\t%s\t%s\n",opt,option_index,optind,optarg,long_options[option_index].name);
	switch(opt){
	    case 'h':
                printf("Usage: test --name your-name --age your-age\n");
                return 0;
                break;
	    case 'n':
		    strcpy(name,optarg);
	        break;
        case 'a':
                age = atoi(optarg);
                break;
        default:
		printf("Unrecognized argument : %c.\n",opt);
  }
}

printf("Your name is %s, your age is %d\n",name,age);

return 0;
}
