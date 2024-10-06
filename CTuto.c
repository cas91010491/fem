#include <stdio.h>
int main(void){
    int **ptr, i, rows = 5, cols =5;
    ptr = (int **) malloc (rows * sizeof(int));

    for (i = 0; i < 5; i++)
        ptr[i] = (int *) malloc(cols * sizeof(int));

    for (int r = 0; r<rows; r++){
        for (int c = 0; c<cols; c++){
            ptr[r][c] = c;
        }
    }
    for (int r = 0; r<rows; r++)
        for (int c = 0; c<cols; c++)
            printf("%d \t", ptr[r][c]);


    int sum = 0;

    for (i=0;i<5;i++)
        sum += ptr[i];

    printf("sum = %d",sum);

    return 0;

}