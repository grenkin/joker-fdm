#ifndef VAR_EXPR_H_INCLUDED
#define VAR_EXPR_H_INCLUDED

struct VarExpr {
    int num;
    std::vector<int> ind;
    std::vector<double> val;
    double rhs;

    VarExpr ()
        : num(0), rhs(0) {}
    VarExpr (double a)
        : num(0), rhs(-a) {}
    VarExpr& operator*= (double k)
    {
        for (int i = 0; i < num; ++i)
            val[i] *= k;
        rhs *= k;
        return *this;
    }
    VarExpr& operator/= (double k)
    {
        for (int i = 0; i < num; ++i)
            val[i] /= k;
        rhs /= k;
        return *this;
    }
    VarExpr& operator+= (const VarExpr& v)
    {
        for (int vi = 0; vi < v.num; ++vi) {
            bool f = false;
            for (int i = 0; i < num; ++i) {
                if (ind[i] == v.ind[vi]) {
                    f = true;
                    val[i] += v.val[vi];
                    break;
                }
            }
            if (!f) {
                ++num;
                ind.resize(num);
                val.resize(num);
                ind[num - 1] = v.ind[vi];
                val[num - 1] = v.val[vi];
            }
        }
        rhs += v.rhs;
        return *this;
    }
    VarExpr& operator-= (const VarExpr& v)
    {
        for (int vi = 0; vi < v.num; ++vi) {
            bool f = false;
            for (int i = 0; i < num; ++i) {
                if (ind[i] == v.ind[vi]) {
                    f = true;
                    val[i] -= v.val[vi];
                    break;
                }
            }
            if (!f) {
                ++num;
                ind.resize(num);
                val.resize(num);
                ind[num - 1] = v.ind[vi];
                val[num - 1] = - v.val[vi];
            }
        }
        rhs -= v.rhs;
        return *this;
    }
};

VarExpr operator* (const VarExpr& v, double k)
{
    VarExpr r = v;
    return r *= k;
}

VarExpr operator* (double k, const VarExpr& v)
{
    return v * k;
}

VarExpr operator/ (const VarExpr& v, double k)
{
    VarExpr r = v;
    return r /= k;
}

VarExpr operator+ (const VarExpr& u, const VarExpr& v)
{
    VarExpr r = u;
    return r += v;
}

VarExpr operator- (const VarExpr& u, const VarExpr& v)
{
    VarExpr r = u;
    return r -= v;
}

#endif // VAR_EXPR_H_INCLUDED
