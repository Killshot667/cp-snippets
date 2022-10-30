// matrix expo code. 3 snippets (the parts)

vector<vector<int>> multiply(vector<vector<int>> &a, vector<vector<int>> &b)
{
    int p = a.size();
    int q = a[0].size();
    int r = b[0].size();
    vector<vector<int>> result(p, vector<int>(r, 0));
    for (int i = 0; i < p; i++)
    {
        for (int j = 0; j < r; j++)
        {
            for (int k = 0; k < q; k++)
            {
                result[i][j] += (a[i][k] * b[k][j]);
            }
        }
    }
    return result;
}
 
vector<vector<int>> multiplyMod(vector<vector<int>> &a, vector<vector<int>> &b)
{
    int p = a.size();
    int q = a[0].size();
    int r = b[0].size();
    vector<vector<int>> result(p, vector<int>(r, 0));
    for (int i = 0; i < p; i++)
    {
        for (int j = 0; j < r; j++)
        {
            for (int k = 0; k < q; k++)
            {
                result[i][j] = (result[i][j] + a[i][k] * b[k][j] % mod) % mod;
            }
        }
    }
    return result;
}
 
vvi expo(vvi x, int y)
{
    int n = x.size();
    vvi p(n, vi(n, 0));
    rep(i, 0, n)
    p[i][i] = 1;
 
    while (y > 0)
    {
        if (y % 2)
        {
            p = multiplyMod(p, x);
        }
        x = multiplyMod(x, x);
        y = y / 2;
    }
    return p;
}
 






// All possible (nCr)%m value in O(logm) using O(n) precomputation

factorial[0] = 1;
for (int i = 1; i <= MAXN; i++) {
    factorial[i] = factorial[i - 1] * i % m;
}

int nCr(int n, int r) {
    return factorial[n] * modInverse(factorial[r] * factorial[n - r] % mod,mod) % mod;
}

int modInverse(int a, int m)
{
    int x, y;
    int g = ext_gcd(a, m, x, y);
    assert(g == 1);
    int res = (x % m + m) % m;
    return res;
}


// euclid's algorithm

int gcd(int a, int b)
{
    if (b == 0)
    {
        return a;
    }
    return gcd(b, a % b);
}

// extended euclid's algorithm

int ext_gcd(int a, int b, int &x, int &y)
{
    if (b == 0)
    {
        x = 1;
        y = 0;
        return a;
    }

    int x1, y1;
    int d = ext_gcd(b, a % b, x1, y1);
    x = y1;
    y = x1 - y1 * (a / b);
    return d;
}

// finding some solution to linear diophantine equation
bool find_any_solution(int a, int b, int c, int &x0, int &y0, int &g) {
    g = ext_gcd(abs(a), abs(b), x0, y0);
    if (c % g) {
        return false;
    }

    x0 *= c / g;
    y0 *= c / g;
    if (a < 0) x0 = -x0;
    if (b < 0) y0 = -y0;
    return true;
}



// modular iverse for all numbers from 1 to m-1

inv[1] = 1;
for (int i = 2; i < m; ++i)
    inv[i] = m - (m / i) * inv[m % i] % m;

// calculate pow(a,b)
int power(int x, int y)
{
    int temp;
    if ( y == 0)
        return 1;
    temp = power(x, y / 2);
    if (y % 2 == 0)
        return temp * temp;
    else
        return x * temp * temp;
}

/* Iterative Function to calculate (x^y) in O(log y) */
int power(int x, int y)
{
    int res = 1;     // Initialize result

    while (y > 0)
    {
        // If y is odd, multiply x with result
        if (y & 1)
            res = res * x;

        // y must be even now
        y = y >> 1; // y = y/2
        x = x * x; // Change x to x^2
    }
    return res;
}

/* Iterative Function to calculate (x^y)%p in O(log y) */
int powerMod(int x, int y, int p)
{
    int res = 1;     // Initialize result

    x = x % p; // Update x if it is more than or
    // equal to p

    if (x == 0) return 0; // In case x is divisible by p;

    while (y > 0)
    {
        // If y is odd, multiply x with result
        if (y & 1)
            res = (res * x) % p;

        // y must be even now
        y = y >> 1; // y = y/2
        x = (x * x) % p;
    }
    return res;
}

// reversing a number
int reversDigits(int num)
{
    int rev_num = 0;
    while (num > 0)
    {
        rev_num = rev_num * 10 + num % 10;
        num = num / 10;
    }
    return rev_num;
}

// checking number is palindrome

bool isPalindrome(int n)
{
    // Find the appropriate divisor
    // to extract the leading digit
    int divisor = 1;
    while (n / divisor >= 10)
        divisor *= 10;

    while (n != 0)
    {
        int leading = n / divisor;
        int trailing = n % 10;

        // If first and last digit
        // not same return false
        if (leading != trailing)
            return false;

        // Removing the leading and trailing
        // digit from number
        n = (n % divisor) / 10;

        // Reducing divisor by a factor
        // of 2 as 2 digits are dropped
        divisor = divisor / 100;
    }
    return true;
}

// cobmination using pascal triangle and dp
int nCrModp(int n, int r, int p)
{
    // Optimization for the cases when r is large
    if (r > n - r)
        r = n - r;

    // The array C is going to store last row of
    // pascal triangle at the end. And last entry
    // of last row is nCr
    int C[r + 1];
    memset(C, 0, sizeof(C));

    C[0] = 1; // Top row of Pascal Triangle

    // One by constructs remaining rows of Pascal
    // Triangle from top to bottom
    for (int i = 1; i <= n; i++) {

        // Fill entries of current row using previous
        // row values
        for (int j = min(i, r); j > 0; j--)

            // nCj = (n-1)Cj + (n-1)C(j-1);
            C[j] = (C[j] + C[j - 1]) % p;
    }
    return C[r];
}

// Combination using Lucas's theorem
// Lucas Theorem based function that returns nCr % p
// This function works like decimal to binary conversion
// recursive function.  First we compute last digits of
// n and r in base p, then recur for remaining digits
int nCrModpLucas(int n, int r, int p)
{
    // Base case
    if (r == 0)
        return 1;

    // Compute last digits of n and r in base p
    int ni = n % p, ri = r % p;

    // Compute result for last digits computed above, and
    // for remaining digits.  Multiply the two results and
    // compute the result of multiplication in modulo p.
    return (nCrModpLucas(n / p, r / p, p) * // Last digits of n and r
            nCrModpDP(ni, ri, p)) % p;  // Remaining digits
}


void sieve(int n)
{

    bool prime[n + 1];
    memset(prime, true, sizeof(prime));
    prime[0] = prime[1] = false;

    for (int p = 2; p * p <= n; p++)
    {

        if (prime[p] == true)
        {

            for (int i = p * p; i <= n; i += p)
                prime[i] = false;
        }
    }

}


// euler totient function in sqrt(n)

int ETF1(int n)
{
    int res = n;
    for (int i = 2; i * i <= n; i++)
    {
        if (n % i == 0)
        {
            res /= i;
            res *= (i - 1);
            while (res % i == 0)
            {
                res /= i;
            }
        }
    }
    if (res > 1)
    {
        res /= n;
        res *= (n - 1);
    }
}

// euler totient for 1 to n in nlog(log(n))
int N = 1000002;
int phi[1000004];

void ETF2()
{
    for (int i = 1; i <= N; i++)
        phi[i] = i;

    for (int i = 2; i <= N; i++)
    {
        if (phi[i] == i)
        {
            for (int j = i; j <= N; j += i)
            {
                phi[j] /= i;
                phi[j] *= (i - 1);
            }

        }
    }

}


// Mobius function values usig linear sieve

mob[1] = 1;
for (int i = 2; i < MAXN; ++i) {
    if (!lp[i]) for (int j = i; j < MAXN; j += i)
            if (!lp[j]) lp[j] = i;
    m[i] = [](int x) {
        int cnt = 0;
        while (x > 1) {
            int k = 0, d = lp[x];
            while (x % d == 0) {
                x /= d;
                ++k;
                if (k > 1) return 0;
            }
            ++cnt;
        }
        if (cnt & 1) return -1;
        return 1;
    }(i);
}