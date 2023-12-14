#include <fstream>
#include <cmath>
#include <iostream>
#include <iomanip>
const double e = 0.000799725808;
const long double Pi = 3.14159265358;
const double eps = 0.00001;
const double a = 4376;
const double mu = 3.9837242 * pow(10, 14);

double rad_Find(double teta, double p)
{
    return p / (1 + (e * cos(teta)));
}

double speed_Rad(double p, double teta)
{
    return sqrt(mu / p) * e * sin(teta);
}

double speed_N(double p, double teta)
{
    return sqrt(mu / p) * (1 + e * cos(teta));
}

double full_Speed(double p, double teta)
{
    double rad = speed_Rad(p, teta);
    double n = speed_N(p, teta);
    return sqrt(rad * rad + n * n);
}

double excentric_To_True(double E)
{
    if (atan(sqrt((1 + e) / (1 - e)) * tan(E / 2)) * 2 > 0)
    {
        return atan(sqrt((1 + e) / (1 - e)) * tan(E / 2)) * 2;
    }
    else
    {
        return atan(sqrt((1 + e) / (1 - e)) * tan(E / 2)) * 2 + 2 * Pi;
    }
}

double iteration_Method(double Enext, double Enow, double M)
{

    if (fabs(Enow - Enext) < eps)
    {

        return Enext;//запись в файл
    }
    else
    {
        return iteration_Method(e * sin(Enext) + M, Enext, M);
    }
}

int main()
{
    std::cout << "Distance to Earth " << "\t" << "Radial speed " << "\t" << "Normal speed " << "\t" << "Full speed " << std::endl;
    for (int i = 0; i < 361; i++)
    {
        double p = a * (1 - e * e);
        double teta = excentric_To_True(iteration_Method(e * sin(i * 2 * Pi / 360) + (i * 2 * Pi / 360), i * 2 * Pi / 360, i * 2 * Pi / 360));
        std::cout << std::fixed << rad_Find(teta, p) << "\t" << speed_Rad(p, teta) << "\t" << speed_N(p, teta) << "\t" << full_Speed(p, teta) << std::endl;
    }

    return 0;
}