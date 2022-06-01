#include <stdio.h>
#include <math.h>

#define Height_Normal(x, y, x1, y1, x2, y2) (((x2-x1)*(y-y1)-(y2-y1)*(x-x1))/sqrt(pow(x2-x1,2)+pow(y2-y1,2)))
#define Eps (1e-5)/2

int main(void) {
    double m, n, r, x, y;
    long long name = 1;
    scanf("%lf %lf %lf %lf %lf", &m, &n, &r, &x, &y);

    double h12, h23, h34, h45, h56, h61;
    double halfHeight = r * sqrt(3) / 2;

    for (long long iy = 1; (iy + 1) * halfHeight <= n; ++iy) {
        for (long long ix = ((iy % 2) ? 2 : 5); (ix + 2) * r / 2 <= m; ix += 6) {
            h12 = Height_Normal(x, y, (ix - 2) * r / 2, (iy) * halfHeight, (ix - 1) * r / 2, (iy + 1) * halfHeight);
            h23 = Height_Normal(x, y, (ix - 1) * r / 2, (iy + 1) * halfHeight, (ix + 1) * r / 2, (iy + 1) * halfHeight);
            h34 = Height_Normal(x, y, (ix + 1) * r / 2, (iy + 1) * halfHeight, (ix + 2) * r / 2, (iy) * halfHeight);
            h45 = Height_Normal(x, y, (ix + 2) * r / 2, (iy) * halfHeight, (ix + 1) * r / 2, (iy - 1) * halfHeight);
            h56 = Height_Normal(x, y, (ix + 1) * r / 2, (iy - 1) * halfHeight, (ix - 1) * r / 2, (iy - 1) * halfHeight);
            h61 = Height_Normal(x, y, (ix - 1) * r / 2, (iy - 1) * halfHeight, (ix - 2) * r / 2, (iy) * halfHeight);
            if ((fabs(h12) <= Eps) || (fabs(h23) <= Eps) || (fabs(h34) <= Eps) ||
                    (fabs(h45) <= Eps) || (fabs(h56) <= Eps) || (fabs(h61) <= Eps)) {
                printf("EDGE");
                return 0;
            } else if ((h12 < -Eps) & (h23 < -Eps) & (h34 < -Eps) &
                    (h45 < -Eps) & (h56 < -Eps) & (h61 < -Eps)) {
                printf("%llu", name);
                return 0;
            }
            ++name;
        }
    }
    printf("NONE");
    return 0;
}