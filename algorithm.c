#include <stdio.h>

// 初期条件関数 f(x)
float initialCondition(float x) {
    return 2 * x * (1 - x);
}

// 境界条件関数 g(x)
float boundaryCondition(float x) {
    return 0.0;
}

int main() {
    // パラメータ設定
    float T = 100.0;      // 時間上限
    int M = 1000;         // 時間刻み数
    float L = 1.0;        // 空間幅
    int N = 100;          // 空間刻み数
    float c = 1.0;        // 波の速さ
    float h = L / N;      // 空間刻み幅
    float tau = T / M;    // 時間刻み幅
    float lamb = c * tau / h; // 安定性条件を満たすための係数

    // 配列の宣言と初期化
    float x[N+2];         // 空間座標
    float u[N+2];         // u(x, t)
    float old_u[N+2];     // 前の時刻のu(x, t)
    float new_u[N+2];     // 次の時刻のu(x, t)

    // 初期値設定
    for (int j = 1; j <= N; j++) {
        x[j] = j * h - 0.5 * h;
        old_u[j] = initialCondition(x[j]);
    }
    old_u[0] = -old_u[1];
    old_u[N+1] = -old_u[N]; //境界条件

    for (int j = 1; j <= N; j++) {
        x[j] = j * h - 0.5 * h;
        u[j] = old_u[j] + tau * boundaryCondition(x[j]) + (lamb * lamb / 2) * (old_u[j-1] - 2 * old_u[j] + old_u[j+1]);
    }

    u[0] = -u[1];
    u[N+1] = -u[N]; //境界条件

    // 計算と動画描画
    for (int k = 0; k < M; k++) {
        // 現在の状態を表示
        for (int j = 1; j <= N; j++) {
            printf("(%f, %f) ", x[j], u[j]);
        }
        printf("\n");

        // 新しい時刻のu(x, t)を計算
        for (int j = 1; j <= N; j++) {
            new_u[j] = 2 * u[j] - old_u[j] + lamb * lamb * (u[j-1] - 2 * u[j] + u[j+1]);
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
        }
        u[0] = -u[1];
        u[N+1] = -u[N];
    }

    return 0;
}
