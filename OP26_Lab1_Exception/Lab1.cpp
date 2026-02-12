#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <iomanip> 

using namespace std;

struct FileNotFoundException
{
    string fName;
    FileNotFoundException(string name) : fName(name) {}
};

struct TriggerAlgorithmTwo {};
struct ResetZero {};
struct ResetOne {};
struct Point2D
{
    double arg;
    double val;
};

double LinearInterp(double targetX, const vector<Point2D>& points)
{
    if (points.empty()) return 0.0;

    for (size_t i = 0; i < points.size() - 1; ++i)
    {
        bool inRange = (targetX >= points[i].arg && targetX <= points[i + 1].arg) ||
            (targetX <= points[i].arg && targetX >= points[i + 1].arg);

        if (inRange)
        {
            double x1 = points[i].arg;
            double y1 = points[i].val;
            double x2 = points[i + 1].arg;
            double y2 = points[i + 1].val;

            if (abs(x2 - x1) < 1e-9) return y1; 
            return y1 + (y2 - y1) * (targetX - x1) / (x2 - x1);
        }
    }
    if (targetX < points.front().arg) return points.front().val;
    return points.back().val;
}

double ReadDataFile(string path, double x)
{
    ifstream ifs(path);
    if (!ifs.is_open())
    {
        throw FileNotFoundException(path);
    }

    vector<Point2D> dataset;
    double currentX, currentY;

    if (ifs.peek() != '-' && !isdigit(ifs.peek()))
    {
        string dummy;
        getline(ifs, dummy);
    }

    while (ifs >> currentX >> currentY)
    {
        dataset.push_back({ currentX, currentY });
    }

    return LinearInterp(x, dataset);
}

double FindTextValue(string key)
{
    ifstream ifs("dat3.dat");
    if (!ifs.is_open()) return 0.0;

    string word;
    double number;

    if (isalpha(ifs.peek()))
    {
        string header;
        getline(ifs, header);
    }

    while (ifs >> word >> number)
    {
        if (word == key) return number;
    }
    return 0.0;
}

double MathY(double x)
{
    double expr = 100.0 - x * x;
    if (expr < 0) throw ResetZero(); 

    double res = sqrt(expr);
    if (x * res < 1.0) throw ResetOne(); 

    return log(x * res);
}

double MathYsm(double x, double y)
{
    return MathY(x) * y + 0.7 * MathY(y);
}

double MathTsm(double x, double y)
{
    double subExpr = 4.0 * pow(y, 4) - x * x;
    if (subExpr < 0) return 0.0;

    double mainExpr = 5.0 * pow(x, 4) - 3.0 * x * x + 2.0 * sqrt(subExpr);
    if (mainExpr <= 0) throw ResetZero();

    return log2(mainExpr) * MathYsm(y, x);
}

double CalculateMts(double x, double y)
{
    return x * MathTsm(x, y) - y * MathTsm(x, x);
}

double CalculateM(double x, double y, double z)
{
    try
    {
        return x * CalculateMts(x, y) + z * CalculateMts(z, y);
    }
    catch (ResetZero&)
    {
        return 0.0;
    }
    catch (ResetOne&)
    {
        return 1.0;
    }
}

double FuncU1(double arg) { return atan(asin(sin(3 * arg))); }
double FuncT1(double arg) { return atan(acos(sin(2 * arg))); }

double CalcQqnl(double x, double y, double z)
{
    return x / FuncU1(x) + y * FuncT1(y) - FuncU1(z) * FuncT1(z);
}

double CalcQnkl(double x, double y)
{
    return 1.15 * CalcQqnl(x, y, x + y) - 0.95 * CalcQqnl(y, x, x - y);
}

double ExecuteAlg2(double x, double y)
{
    double p1 = CalcQnkl(x, y);
    double p2 = CalcQnkl(y, x);
    return x * p1 + y * p2 - 0.05 * p1 * p2;
}

double GetU(double x) { return ReadDataFile("dat1.dat", x); }
double GetT(double x) { return ReadDataFile("dat2.dat", x); }

double CalcQkn(double x, double y)
{
    return x / GetU(x) + y * GetT(y);
}

double ProcessQnk(double x, double y, double z)
{
    if (abs(x) <= 5.0) throw TriggerAlgorithmTwo();
    if (abs(x) <= 10.0) throw TriggerAlgorithmTwo();

    try
    {
        return CalcQkn(x, y) + x * CalcQkn(y, z);
    }
    catch (FileNotFoundException&)
    {
        throw TriggerAlgorithmTwo();
    }
}

double SafeQnk(double x, double y, double z)
{
    try
    {
        return ProcessQnk(x, y, z);
    }
    catch (TriggerAlgorithmTwo&)
    {
        return ExecuteAlg2(x, y);
    }
}

double CalcRsv(double x, double y, double z)
{
    if (z > x && z > y)      return z * SafeQnk(x, y, z) - x * y;
    else if (x > y && x > z) return x * SafeQnk(z, y, x) + y * z;
    else if (y > x && y > z) return y * SafeQnk(x, z, y) + x * z;

    return z * SafeQnk(y, z, x) - SafeQnk(z, x, y);
}


double GetGText(string txt) { return FindTextValue(txt); }

double CalcStext(double val, string txt)
{
    if (val > 0)
    {
        return FindTextValue(txt) + val;
    }
    else if (txt.empty())
    {
        return FindTextValue("tet") + GetGText("set") - val;
    }
    else
    {
        return FindTextValue("get") + GetGText(txt);
    }
}

double FindMin(double n1, double n2, double n3, double n4)
{
    return min({ n1, n2, n3, n4 });
}

double CalculateK(double x, double y, double z, string txt)
{
    double minVal;
    if (z > 0) minVal = FindMin(x, y, x - z, y - z);
    else       minVal = FindMin(x, y, z - x, z - y);

    return CalcStext(minVal, txt);
}

double FinalFunction(double r, double m, double k)
{
    return 10.0 * k * pow(r, 2) - m * r;
}

int main()
{
    double varX, varY, varZ;
    string userText;

    cout << "Input X, Y, Z: ";
    if (!(cin >> varX >> varY >> varZ))
    {
        cerr << "Error: Incorrect number format." << endl;
        return 1;
    }

    cout << "Input text string: ";
    cin >> userText;

    double valR = CalcRsv(varX, varY, varZ) +
        0.5 * CalcRsv(varY, varZ, varX) * CalcRsv(varZ, varX, varY);

    double valM = CalculateM(varX, varY, varZ);
    double valK = CalculateK(varX, varY, varZ, userText);

    double result = FinalFunction(valR, valM, valK);

    cout << "\n-----------------------" << endl;
    cout << fixed << setprecision(4); 
    cout << "Computed R: " << valR << endl;
    cout << "Computed M: " << valM << endl;
    cout << "Computed K: " << valK << endl;
    cout << "Function Result: " << result << endl;
    cout << "-----------------------" << endl;

    return 0;
}