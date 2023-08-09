#include <stdio.h>

int f(float x)
{
    int value;
    value = 2 * x * (1 - x);
    return value;
}

int main(void)
{
    int T, M, L, N, c;
    int k, j;
    float h, lamb, tau;
    
    L = 1;
    c = 1;
    N = 100;
    M = 150;
    tau = 1;
    
    h = (float)L / N;
    lamb = tau * c / h;

    printf("λ = %f\n", lamb);

    float u[N+2];
    float old_u[N+2];
    float new_u[N+2];

    /* 初期値設定 */
    for (j = 1; j <= N; j++) {
        float xj;
        xj = j * h - 0.5 * h;
        old_u[j] = f(xj);
    }

    old_u[0] = -old_u[1];
    old_u[N+1] = -old_u[N]; /* 境界条件 */

    for (k = 0; k < M; k++) {
        for (j = 1; j <= N; j++) {
            new_u[j] = 2 * old_u[j] - u[j] + lamb * lamb * (old_u[j-1] - 2 * old_u[j] + old_u[j+1]);
        }

        for (j = 1; j <= N; j++) {
            old_u[j] = u[j];
        }

        old_u[0] = -old_u[1];
        old_u[N+1] = -old_u[N]; /* 境界条件 */

        for (j = 1; j <= N; j++) {
            u[j] = new_u[j];
        }

        u[0] = -u[1];
        u[N+1] = -u[N]; /* 境界条件 */
    }

    return 0;
}
