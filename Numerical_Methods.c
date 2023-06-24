#include <stdio.h>
#include <math.h>
#define MAX 20

float f(float pol[MAX], int p , float y);
float fturev(float pol[MAX], int p , float y);
void Bisection(float a, float b, float pol[], float eps, int p);
void RegulaFalsi(float a, float b, float pol[], float eps, int p);
void NewtonRaphson(float x0, float pol[], float eps, int p);
float determinantal(float [][25], float);
void kofaktor(float [][25], float);
void transpozeal(float [][25], float [][25], float);
void GaussEleminasyon(float A[][20] , int n );
void GaussSeidel(float A[][MAX], float x, float y, float z, float eps);
void SayisalTurev(float pol[MAX], int p , float x , float h);
void Trapez(float pol[MAX], int p , float x0 , float xn , float n);
void Simpson(float pol[MAX], int p , float x0 , float xn , float n);
void gregorynewton(float F[][MAX] , float n , float k);

int main(){
	int w , u , i , j , p , m;
	float eps;
	float x0 , xn , n , l , d , a , b , r , t , y , z , k , x , h;
	
	do{
		float pol[10]={0};
		float B[25][25]={0};
		float G[20][20]={0};
		float A[MAX][MAX]={0};
		float F[MAX][MAX]={0};
		
		printf("Menu :\n");
		printf("(1) Dogrusal Olmayan Esitliklerin Cozumu \n(2) Matris Islemleri \n(3) Sayisal Turev \n(4) Sayisal Integral \n(5) Gregory Newton Enterpolasyonu \n(6) Cikis \n");
		scanf("%d" , &w);
		
		if(w != 6){
			if(w == 1){
				printf("Dogrusal Olmayan Esitliklerin Cozumu \nAlt Menu : \n");
				printf("(1) Bisection \n(2) Regula Falsi \n(3) Newton Raphson \n");
				scanf("%d" , &u);
				if(u == 1){
					printf("Bisection\n\n");
					
					                                                                                                       //1.Yazýlým
					printf("fonksiyonun en buyuk derecesini giriniz : ");
					scanf("%d" , &p);
	
					for(i=p ; i>=0 ; i--){
						printf("%d. dereceden elemanin katsayisini yaziniz : " , i);
						scanf("%f" , &pol[i]);
					}
	
					do{
						printf("Arasinda kok olan aralik seciniz a-b \n");
						printf("a : ");
						scanf("%f" , &a);
		
						printf("b : ");
						scanf("%f" , &b);
		
					}while( (f(pol , p , a) * f(pol , p , b) ) > 0);
	
					printf("eps hata degeri giriniz : ");
					scanf("%f" , &eps);
	
	                Bisection( a , b , pol , eps , p);                                                                                                      
					
					
				}else if(u == 2){
					printf("Regula Falsi\n\n");
					
					                                                                                                        //2.Yazýlým
					printf("fonksiyonun en buyuk derecesini giriniz : ");
					scanf("%d" , &p);
	
					for(i=p ; i>=0 ; i--){
						printf("%d. dereceden elemanin katsayisini yaziniz : " , i);
						scanf("%f" , &pol[i]);
					}
	
					do{
						printf("Arasinda kok olan aralik seciniz a-b \n");
						printf("a : ");
						scanf("%f" , &r);
		
						printf("b : ");
						scanf("%f" , &t);
		
					}while( (f(pol , p , r) * f(pol , p , t) ) > 0);
	
					printf("eps hata degeri giriniz : ");
					scanf("%f" , &eps);
	
					RegulaFalsi( r , t , pol , eps , p);                                                                                                        
					
					
				}else{
					printf("Newton Raphson\n\n");
					
					                                                                                                        //3.Yazýlým
					                                                                                                        
					printf("fonksiyonun en buyuk derecesini giriniz : ");
					scanf("%d" , &p);
	
					for(i=p ; i>=0 ; i--){
						printf("%d. dereceden elemanin katsayisini yaziniz : " , i);
						scanf("%f" , &pol[i]);
					}
	
					printf("Baslangic indisi giriniz x0 : ");
					scanf("%f" , &x0);
	
					printf("eps hata degeri giriniz : ");
					scanf("%f" , &eps);
	
					NewtonRaphson( x0 , pol , eps , p);
					
					
				}
			}else if(w == 2){
				printf("Matris Islemleri \n Alt Menu : \n");
				printf("(1) Matrisin Tersini Alma \n(2) Gauss Eleminasyon \n(3) Gauss Seidel \n");
				scanf("%d" , &u);
				if(u == 1){
					printf("Matrisin Tersini Alma\n\n");
					
					                                                                                                         //4.Yazýlým
					printf("Matrisin satir ve sutun sayisini giriniz : ");
					scanf("%f", &l);
  
  					for (i = 0;i < l; i++){
    					for (j = 0;j < l; j++){
    						printf("Satir[%d] , Sutun[%d] degerini gir : " , i , j);
        					scanf("%f", &B[i][j]);
        				}
    				}
					d = determinantal(B, l);
					printf("determinant = %f\n" , d);
					if (d == 0){
						printf("\nMatrisin tersi alinamaz\n");
					}
					else{
  						kofaktor(B, l);
					}                                                                                                         
					                                                                                                         
					                                                                                                         
					
					
				}else if(u == 2){
					printf("Gauss Eleminasyon\n\n");
					
					                                                                                                         //5.Yazýlým
					                                                                                                         
					printf("Matrisin satir ve sutun sayisini giriniz : ");
    				scanf("%d",&m);
    
    				for(i=1; i<=m; i++){
        				for(j=1; j<=(m+1); j++){
            				printf("Satir[%d] , Sutun[%d] degerini gir : " , i , j);
            				scanf("%f",&G[i][j]);
        				}
    				}
    				GaussEleminasyon(G , m);                                                                                                         
					
					
				}else{
					printf("Gauss Seidel\n\n");
					
					                                                                                                          //6.Yazýlým
					printf("3 bilinmeyeni olan denklem sisteminin degerlerini giriniz\n");
	
					for(i=1 ; i<=3 ; i++){
						for(j=1 ; j<=4 ; j++){
							printf("A[%d][%d] : " , i , j);
							scanf("%f" , &A[i][j]);
						}
					}
					printf("x baslangic degeri giriniz : ");
					scanf("%f" , &x);
	
					printf("y baslangic degeri giriniz : ");
					scanf("%f" , &y);
	
					printf("z baslangic degeri giriniz : ");
					scanf("%f" , &z);
	
					printf("epsilon hata degeri giriniz : ");
					scanf("%f" , &eps);
		
					GaussSeidel(A , x , y , z , eps);                                                                                                          
					                                                                                                          
					                                                                                                          
					
					
					
				}
			}else if(w == 3){
				printf("Sayisal Turev \n\n");
				
				                                                                                                               //7.Yazýlým
					printf("fonksiyonun en buyuk derecesini giriniz : ");
					scanf("%d" , &p);
	
					for(i=p ; i>=0 ; i--){
						printf("%d. dereceden elemanin katsayisini yaziniz : " , i);
						scanf("%f" , &pol[i]);
					}
	
					printf("Turevi alinacak x noktasi degeri giriniz : ");
					scanf("%f" , &x);
	
					printf("Turev icin h degeri giriniz : ");
					scanf("%f" , &h);
	
					SayisalTurev(pol , p , x , h);                                                                                                               
				                                                                                                               
				                                                                                                               
				
			}else if(w == 4){
				printf("Sayisal Integral \n Alt Menu : \n");
				printf("(1) Trapez \n(2) Simpson \n");
				scanf("%d" , &u);
				if(u == 1){
					printf("Trapez \n\n");
					
					                                                                                                           //8.Yazýlým
					printf("fonksiyonun en buyuk derecesini giriniz : ");
					scanf("%d" , &p);
		
					for(i=p ; i>=0 ; i--){
						printf("%d. dereceden elemanin katsayisini yaziniz : " , i);
						scanf("%f" , &pol[i]);
					}
					printf("Integrali alinacak araligin x0 ve xn uc degerlerini giriniz\n");
					printf("x0 : ");
					scanf("%f" , &x0);
	
					printf("xn : ");
					scanf("%f" , &xn);
		
					printf("Esit parca sayisini giriniz n : ");
					scanf("%f" , &n);
	
					Trapez(pol , p , x0 , xn , n);                                                                                                           
					                                                                                                           
					                                                                                                           
					
					
				}else{
					printf("Simpson \n\n");
					
					                                                                                                           //9.Yazýlým
					printf("fonksiyonun en buyuk derecesini giriniz : ");
					scanf("%d" , &p);
	
					for(i=p ; i>=0 ; i--){
						printf("%d. dereceden elemanin katsayisini yaziniz : " , i);
						scanf("%f" , &pol[i]);
					}
					printf("Integrali alinacak araligin x0 ve xn uc degerlerini giriniz\n");
					printf("x0 : ");
					scanf("%f" , &x0);
	
					printf("xn : ");
					scanf("%f" , &xn);
	
					printf("Esit parca sayisini giriniz (cift olmali)    n : ");
					scanf("%f" , &n);
		
					Simpson(pol , p , x0 , xn , n);                                                                                                           
					                                                                                                           
					                                                                                                           
					
					
				}
			}else if(w == 5){
				printf("Gregory Newton Enterpolasyonu \n\n");
				
				                                                                                                              //10.Yazýlým
					printf("Enterpolasyon icin x0 - xn degerleri giriniz\n");
					printf("x0 : ");
					scanf("%f" , &x0);
	
					printf("xn : ");
					scanf("%f" , &xn);
	
					printf("h degeri giriniz : ");
					scanf("%f" , &h);
	
					n = (xn - x0) / h;
	
					for(i=0 ; i<n+1 ; i++){
						printf("f(x%d) degerini giriniz : " , i);
						scanf("%f" , &F[i][0]);
					}	
	
					printf("x degeri giriniz : ");
					scanf("%f" , &x);
	
					k = (x - x0) / h;
	
					gregorynewton(F , n , k);                                                                                                              
				                                                                                                              
				                                                                                                              
				
			}
		}
	}while(w != 6);
	
	
return 0;
}

float f(float pol[MAX], int p , float y){
	int i,j;
	float total=0 , c;
	
	for(i=0 ; i<=p ; i++){
		c=pol[i];
		total = total + (c * pow(y,i));
		
	}
	printf("f(%f) : %f\n" , y , total);
	return total;
}

float fturev(float pol[MAX], int p , float y){
	int i,j;
	float total=0 , c;
	
	for(i=1 ; i<=p ; i++){
		c=pol[i];
		total = total + (c * i * pow(y,i-1));
		
	}
	return total;
}

void Bisection(float a, float b, float pol[], float eps, int p){
	int i;
	float high , low , mid;
	
	
	low = a;
	high = b;
	mid = ( low + high ) / 2;
	i=1;
	
	printf("mid : %f\n" , mid);
	
	while( ( ( high - low) / pow(2,i) ) > eps ){
		if( f(pol , p , low) * f(pol , p , mid) < 0 ){
			high = mid;
		}
		else{
			low = mid;
		}
		mid = ( low + high ) / 2;
		i++;
		
		printf("mid : %f\n" , mid);
	}
	printf("kok : %f\n" , mid);
	printf("\n\n");
}

void RegulaFalsi(float a, float b, float pol[], float eps, int p){
	int i;
	float high , low , mid;
	
	low = a;
	high = b;
	mid = ( high*f(pol , p , low) - low*f(pol , p , high) ) / ( f(pol , p , low) - f(pol , p , high) );
	i = 1;
	
	printf("mid : %f\n" , mid);
	
	while( ( ( high - low) / pow(2,i) ) > eps ){
		if( f(pol , p , low) * f(pol , p , mid) < 0 ){
			high = mid;
		}
		else{
			low = mid;
		}
		mid = ( high*f(pol , p , low) - low*f(pol , p , high) ) / ( f(pol , p , low) - f(pol , p , high) );
		i++;
		
		printf("mid : %f\n" , mid);
	}
	
	printf("kok = %f\n" , mid);
	printf("\n\n");
	
}

void NewtonRaphson(float x0, float pol[], float eps, int p){
	float x1;
	
	x1 = x0 - ( f(pol , p , x0) / fturev(pol , p , x0) );
	
	printf("x1 : %f\n" , x1);
	
	while( fabs(x1 - x0) > eps){
		x0 = x1;
		x1 = x0 - ( f(pol , p , x0) / fturev(pol , p , x0) );
		printf("x1 : %f\n" , x1);
	}
	
	printf("kok : %f\n" , x1);
	printf("\n\n");
}

float determinantal(float a[25][25], float k){
	
	float s = 1, det = 0, b[25][25];
	int i, j, m, n, c;
	
	if (k == 1){
    	return (a[0][0]);
    }
	else{
    	det = 0;
    	for (c = 0; c < k; c++){
        	m = 0;
        	n = 0;
        	for (i = 0;i < k; i++){
            	for (j = 0 ;j < k; j++){
            		b[i][j] = 0;
                	if (i != 0 && j != c){
                		b[m][n] = a[i][j];
                		if (n < (k - 2)){
                			n++;
						}
                    	else{
                    		n = 0;
                    		m++;
                		}
        			}
        		}
        	}
        	det = det + s * ( a[0][c] * determinantal(b, k - 1) );
        	s = -1 * s;
		}
    }
 
return (det);
}

void kofaktor(float num[25][25], float f){
	
	float b[25][25], fac[25][25];
	int p, q, m, n, i, j;
 
	for (q = 0;q < f; q++){
		for (p = 0;p < f; p++){
    		m = 0;
    		n = 0;
    		for (i = 0;i < f; i++){
    			for (j = 0;j < f; j++){
        			if (i != q && j != p){
            			b[m][n] = num[i][j];
            			if (n < (f - 2)){
            				n++;
            			}
            			else{
               				n = 0;
               				m++;
            			}
        			}
        		}
    		}
    		fac[q][p] = pow(-1, q + p) * determinantal(b, f - 1);
    	}
	}
	transpozeal(num, fac, f);
}

void transpozeal(float num[25][25], float fac[25][25], float r){
	
	int i, j;
	float b[25][25], ters[25][25], d;
 
	for (i = 0;i < r; i++){
    	for (j = 0;j < r; j++){
			b[i][j] = fac[j][i];
        }
    }
    
    printf("Adjoint Matris :\n");
    
    for (i = 0;i < r; i++){
		for (j = 0;j < r; j++){
			printf("\t%f", b[i][j]);
		}
		printf("\n");
	}
    
	d = determinantal(num, r);
	for (i = 0;i < r; i++){
    	for (j = 0;j < r; j++){
    		ters[i][j] = b[i][j] / d;
        }
    }
	printf("\n\n\nMatrisin Tersi : \n");
 
	for (i = 0;i < r; i++){
    	for (j = 0;j < r; j++){
    		printf("\t%f", ters[i][j]);
        }
    	printf("\n");
    }
    printf("\n\n");
    
}

void GaussEleminasyon(float A[][20] , int n ){
	
	int i , j , k;
	float c , toplam = 0.0 , x[10];
	
	for(j=1; j<=n; j++){ /* Alt Ucgen tarafýný sýfýrlama */
    	for(i=1; i<=n; i++){
        	if(i>j){
            	c=A[i][j]/A[j][j];
                for(k=1; k<=n+1; k++){
                	A[i][k]=A[i][k]-c*A[j][k];
                }
            }
        }
    }
    
    for (i = 1;i <= n; i++){
		for (j = 1;j <= n+1; j++){
			printf("\t%f", A[i][j]);
		}
		printf("\n");
	}
    
    x[n]=A[n][n+1]/A[n][n];   /* Son x'in degerini yazma */
    
    for(i=n-1; i>=1; i--){     /* diger x degerlerini bulma */
    	toplam=0;
        for(j=i+1; j<=n; j++){
        	toplam=toplam+A[i][j]*x[j];
        }
    	x[i]=(A[i][n+1]-toplam)/A[i][i];
    }
    
    printf("\nCozum : \n");
    
    for(i=1; i<=n; i++){
    	printf("\nx%d=%f\t",i,x[i]);
    }
    printf("\n\n");
}

void GaussSeidel(float A[][MAX], float x, float y, float z, float eps){
	float dx , dy , dz;
	float a , b , c;
	
	dx=100;
	dy=100;
	dz=100;
	
	while( (dx>eps) || (dy>eps) || (dz>eps) ){
		
		a = (A[1][4] - A[1][2]*y - A[1][3]*z ) / A[1][1];
		dx = fabs(a - x);
		x = a;
		
		b = (A[2][4] - A[2][1]*x - A[2][3]*z ) / A[2][2];
		dy = fabs(b - y);
		y = b;
		
		c = (A[3][4] - A[3][1]*x - A[3][2]*y ) / A[3][3];
		dz = fabs(c - z);
		z = c;
		
		printf("x = %f , y = %f , z = %f  dx = %f , dy = %f , dz = %f \n" , x , y , z , dx , dy , dz);
	}
	printf("\n\n");
}

void SayisalTurev(float pol[MAX], int p , float x , float h){
	float k,l,m;
	
	printf("\n");
	
	k = ( f(pol,p,x+h) - f(pol,p,x) ) / h;
	l = ( f(pol,p,x) - f(pol,p,x-h) ) / h;
	m = ( f(pol,p,x+h) - f(pol,p,x-h) ) / (2*h);
	
	printf("ileri fark turevi : %f\n" , k);
	printf("geri fark turevi : %f\n" , l);
	printf("merkezi fark turevi : %f\n" , m);
	
	printf("\n\n");
}

void Trapez(float pol[MAX], int p , float x0 , float xn , float n){
	int i;
	float h , s;
	
	h = ( xn - x0 ) / n;
	s = 0;
	
	printf("h : %f\n" , h);
	
	for(i=1 ; i<n ; i++){
		s += h * f(pol , p , x0+h*i) ;
	}
	
	s += h * ( ( f(pol , p , x0) + f(pol , p , xn) ) / 2 );
	
	printf("Integral : %f\n" , s);
	
	printf("\n\n");
}

void Simpson(float pol[MAX], int p , float x0 , float xn , float n){
	float h , b , s;
	int i;
	
	h = ( xn - x0 ) / n;
	s = 0;
	
	for(i=1 ; i<n ; i+=2){
		b=h*i;
		s = s + 4 * f(pol , p , x0+b);
	}
	
	for(i=2 ; i<(n-1) ; i+=2){
		b=h*i;
		s = s + 2 * f(pol , p , x0+b);
	}
	
	s = s + f(pol , p , x0) + f(pol , p , xn);
	s = s*(h/3);
	
	printf("Integral : %f\n" , s);
	
	printf("\n\n");
	
}

void gregorynewton(float F[][MAX] , float n ,float k){
	int i , j;
	float z , g , toplam;
	
	for(i=1 ; i<=n ; i++){
		for(j=i ; j<=n ; j++){
			F[j][i] = F[j][i-1] - F[j-1][i-1];
		}
	}
	z = k;                                                     //   z = k.(k-1).(k-2)...
	g = 1;                                                     //   g = 1.2.3.4....
	
	toplam = F[0][0];
	for(i=1 ; i<=n ; i++){
		
		toplam += (z * F[i][i]) / g;
		
		z = z * (k-1);
		k = k - 1;
		g = g * (i+1);
	}
	printf("sonuc : %f\n" , toplam);
	
	printf("\n\n");
		
}



