#include <stdio.h>

// 初期条件関数 f(x)
double initialCondition(double x) {
    return 2.0 * x * (1.0 - x);
}

// 境界条件関数 g(x)
double boundaryCondition(double x) {
    return 0.0;
}

int main() {
    double T = 10.0; // 時間上限
    int M = 10000; // 時間刻み数
    double L = 1.0; // 空間幅(x方向)
    int N = 100; // 空間刻み数(x方向)
    double H = 0.02; // 空間幅(y方向)
    int S = 2; // 空間刻み数(y方向)
    double tau = T / (M * 1.0); // 時間刻み幅
    double alpha = (L * 1.0) / (N * 1.0); // 空間刻み幅(x方向)
    double beta = (H * 1.0) / (S * 1.0); // 空間刻み幅(y方向)

	// 配列の宣言と初期化
    double x[N+2]; // 空間座標(x方向)
    // double y[N+2]; //空間座標(y方向)
    double u[N+2]; // u(x, t)
    double old_u[N+2]; // 前の時刻のu(x, t)
    double new_u[N+2]; // 次の時刻のu(x, t)
    double dif_u[N+2]; // 差
    double plus_u[N+2]; // y軸手前(負)方向
    double minus_u[N+2]; //y軸奥(正)方向
    double t[M];    // 時間


    // 初期値設定
    for (int j = 1; j <= N; j++) {
            x[j] = j * alpha - 0.5 * alpha;
            old_u[j] = initialCondition(x[j]);
        }

	//境界条件
	u[0] = -u[1];
    u[N+1] = -u[N]; //境界条件

    // 計算と動画描画
    for (int k = 0; k < 6; k++) {
        t[k] = ((k + 1) * 1.0) * tau;
        //初期条件
        u[1] = old_u[N-1];
        u[2] = old_u[N];

        //横一列(y方向の影響考えない)ver(★いい感じに動いたら変えるところ★)
        minus_u[1] = u[1];
        plus_u[2] = u[2];

        for (int i = 1; i <= S; i++) {
            // y[i] = i * beta - 0.5 * beta;
            // 現在の状態を表示
            for (int j = 1; j <= N; j++) {
                if (k == 5)
                    printf("%f %f %f \n", x[j], t[k], u[j]);
            }
            if (k == 5)
                printf("\n");

            // 新しい時刻のu(x, t)を計算
            for (int j = 3; j <= N; j++) {
                            //jの前後の波高（u）の差 (new_u[j-1] - new_u[j+1])
                            //u[i+1]→plus_u[j]、u[i-1]→minus_u[j
                            dif_u[j] = old_u[j-1] - old_u[j+1] + (tau / alpha) * ((3.0 / 2.0) * (u[j+1] - u[j-1]) * (u[j+1] - u[j-1]) + 6.0 * u[j] * (u[j+1] - 2.0 * u[j] + u[j-1]) + (1 / alpha * alpha) *(u[j+2] - 4.0 * u[j+1] + 6*u[j] - 4.0 * u[j-1] + u[j-2])) + 3.0 * (alpha * tau /(beta * beta)) * (plus_u[j] - 2.0 * u[j] + minus_u[j]);
                            //(-1) * dif_uを代入
                            new_u[j] = (u[j-1] + (-1.0) * dif_u[j-1]) / 2.0;
            }

            // 前の時刻のu(x, t)を更新
            for (int j = 1; j <= N; j++) {
                old_u[j] = u[j];
            }
            old_u[0] = -old_u[1];
            old_u[N+1] = -old_u[N];

            // u(x, t)を更新
            for (int j = 1; j <= N; j++) {
                u[j] = new_u[j];
                //横一列(y方向の影響考えない)ver(★いい感じに動いたら変えるところ★)
                minus_u[j] = u[j];
                plus_u[j] = u[j];
            }
            u[0] = -u[1];
            u[N+1] = -u[N];
        }
    }
    return 0;
}