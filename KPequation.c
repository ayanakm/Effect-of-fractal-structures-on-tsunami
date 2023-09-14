#include <stdio.h>
#include <math.h>

// 初期条件関数 f(x)
double initialCondition(double x) {
    return expl(-100 * (x - 0.4) * (x - 0.4));
}

// 境界条件関数 g(x)
double boundaryCondition(double x) {
    return 0.0;
}

int main() {
    double T = 0.0; // 時間上限
    int M = 0; // 時間刻み数
    double L = 0.0; // 空間幅(x方向)
    int N = 0; // 空間刻み数(x方向)
    double H = 0.0; // 空間幅(y方向)
    int S = 0; // 空間刻み数(y方向)

    // 初期値設定
    T = 10.0; // 時間上限
    M = 100000; // 時間刻み数
    L = 1.0; // 空間幅(x方向)
    N = 100; // 空間刻み数(x方向)
    H = 0.02; // 空間幅(y方向)
    S = 2; // 空間刻み数(y方向)
    

    double tau = T / (M * 1.0); // 時間刻み幅
    double alpha = (L * 1.0) / (N * 1.0); // 空間刻み幅(x方向)
    double beta = (H * 1.0) / (S * 1.0); // 空間刻み幅(y方向)

	// 配列の宣言と初期化
	double x[N+2]; // 空間座標(x方向)
	double y[S+2]; // 空間座標(y方向)
    double u[S+2][N+2]; // u(x, t)
    double old_u[S+2][N+2]; // 前の時刻のu(x, t)
    double new_u[S+2][N+2]; // 次の時刻のu(x, t)
    double dif_u[S+2][N+2]; // 差
	double t[M]; // 時間

    // 初期条件
    for (int i = 0; i <= S+1; i++){
        for (int j = 0; j <= N+1; j++) {
            x[j] = j * alpha - 0.5 * alpha;
            y[i] = i * beta - 0.5 * beta;
            u[i][j] = initialCondition(x[j]);
            // printf("%f %f %f \n", x[j], y[i], old_u[i][j]);
            u[i][j] = old_u[i][j] + 3;
        }
    }

    for (int i = 0; i <= S+1; i++){
        for (int j = 0; j <= N+1; j++) {
            if (i == 0){
                u[0][j] = -u[1][j];
            }
            if (i == S+1){
                u[S+1][j] = -u[S][j];
            }
            if (j == 0){
                u[i][0] = u[i][1];
            }
            if (j == N+1){
                u[i][N+1] = u[i][N];
            }
            // printf("i = %d, j = %d, u = %10.8e \n", i, j, u[i][j]);
        }
        // printf("\n");
    }

    //二重for文でold_uにuを代入
    for (int i = 0; i <= S+1; i++){
        for (int j = 0; j <= N+1; j++) {
            old_u[i][j] = u[i][j];
        }
    }

    // 計算と動画描画
    for (int k = 0; k < 11; k++) {
        t[k] = ((k + 1) * 1.0) * tau;

        // printf("t = %d \n", k);
        for (int i = 0; i <= S+1; i++){
            // 現在の状態を表示
            for (int j = 1; j <= N; j++) {
                // if (k == 5)
                // printf("i = %d, j = %d, %f %f %10.8e \n",i ,j, x[j], y[i], u[i][j]);
                if (i > 0 && i < S + 1){
                    printf("%f %f %10.8e \n",x[j], y[i], u[i][j]);
                }
            }

            // 新しい時刻のu(x, t)を計算
            for (int j = 3; j <= N; j++) {
                // printf("dif_u = " "%f \n", dif_u[i][j]);
                //jの前後の波高（u）の差 (new_u[j-1] - new_u[j+1])
                dif_u[i][j] = old_u[i][j-1] - old_u[i][j+1] + (tau / alpha) * ((3.0 / 2.0) * (u[i][j+1] - u[i][j-1]) * (u[i][j+1] - u[i][j-1]) + 6.0 * u[i][j] * (u[i][j+1] - 2.0 * u[i][j] + u[i][j-1])) + (1.0 / alpha * alpha) * (u[i][j+2] - 4.0 * u[i][j+1] + 6.0 * u[i][j] - 4.0 * u[i][j-1] + u[i][j-2]) + 3.0 * (alpha * tau /(beta * beta)) * (u[i+1][j] - 2.0 * u[i][j] + u[i-1][j]);
                //(-1) * dif_uを代入
                new_u[i][j] = u[i][j-1] + ((-1.0) * dif_u[i][j-1])/2;
            }

            // 前の時刻のu(x, t)を更新
            for (int j = 1; j <= N; j++) {
                old_u[i][j] = u[i][j];
            }
            old_u[i][0] = -old_u[i][1];
            old_u[i][N+1] = -old_u[i][N];

            // u(x, t)を更新
            for (int j = 3; j <= N; j++) {
                //横一列(y方向の影響考えない)ver(★いい感じに動いたら変えるところ★)
                u[i][j] = new_u[i][j];
            }
            //j=1,2は前の時刻のj=N,N-1の値
            u[i][1] = old_u[i][N-1];
            u[i][2] = old_u[i][N];

            //jの境界条件をつける（最後もう一回2重for文回す）
            //境界条件
            u[i][0] = -u[i][1];
            u[i][N+1] = -u[i][N];


        }

        printf("\n");

        for (int i = 0; i <= S+1; i++){
            for (int j = 0; j <= N+1; j++) {
                if (i == 0){
                    u[0][j] = -u[1][j];
                }
                if (i == S+1){
                    u[S+1][j] = -u[S][j];
                }
                if (j == 0){
                    u[i][0] = u[i][1];
                }
                if (j == N+1){
                    u[i][N+1] = u[i][N];
                }
            }
            // printf("\n");
        }

    }
    return 0;
}  
