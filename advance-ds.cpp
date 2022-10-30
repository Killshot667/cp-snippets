
// default fenwick tree for cumulative sum of int
struct FENWICK {
    int N;
    vector<int> tree;


    // constructor
    void fenwick(int n)
    {
        N = n;
        tree.clear();
        tree.resize(N + 1);
    }

    // updating the element of the array at index idx by val
    void update(int idx, int val)
    {
        while (idx <= N)
        {
            tree[idx] += val;
            idx += (idx & -idx);
        }
    }

    // reading cumulative freq at idx (prefix sum)
    int read(int idx)
    {
        int sum = 0;
        while (idx > 0) {
            sum += tree[idx];
            idx -= (idx & -idx);
        }
        return sum;
    }


    // building fenwick tree. Call this first
    void build(vector<int> a)
    {
        for (int i = 1; i <= n; i++)
            update(i, a[i]);
    }

    // scaling all elements by c (multiply by c)
    void scale(int c)
    {
        for (int i = 1; i <= N; i++)
            tree[i] = tree[i] * c;
    }


    // find the value of a particular array element
    int readSingle(int idx)
    {
        int sum = tree[idx]; // this sum will be decreased
        if (idx > 0)
        {   // the special case
            int z = idx - (idx & -idx);
            idx--; // idx is not important anymore, so instead y, you can use idx
            while (idx != z) { // at some iteration idx (y) will become z
                sum -= tree[idx];
                // substruct tree frequency which is between y and "the same path"
                idx -= (idx & -idx);
            }
        }
        return sum;
    }


    //finding the index with a particular cumulative freq when the entire cumulative freq array is non-decreasing

    // If in the tree exists more than one index with the same
    // cumulative frequency, this procedure will return
    // some random one of them

    // bitMask - initialy, it is the greatest bit of MaxIdx
    // bitMask stores the current interval that should be searched
    int findR(int cumFre)
    {
        int idx = 0; // this variable will be the output
        bitmask = 0;
        int temp = N / 2;
        while (temp != 0)
        {
            temp = temp / 2;
            bitmask++;
        }

        while (bitMask != 0) {
            int tIdx = idx + bitMask; // the midpoint of the current interval
            bitMask >>= 1; // halve the current interval
            if (tIdx > N) // avoid overflow
                continue;
            if (cumFre == tree[tIdx]) // if it is equal, simply return tIdx
                return tIdx;
            else if (cumFre > tree[tIdx]) {
                // if the tree frequency "can fit" into cumFre,
                // then include it
                idx = tIdx; // update index
                cumFre -= tree[tIdx]; // update the frequency for the next iteration
            }
        }
        if (cumFre != 0) // maybe the given cumulative frequency doesn't exist
            return -1;
        else
            return idx;
    }

    // If in the tree exists more than one index with a same
    // cumulative frequency, this procedure will return
    // the greatest one
    int findG(int cumFre) {
        int idx = 0;

        bitmask = 0;
        int temp = N / 2;
        while (temp != 0)
        {
            temp = temp / 2;
            bitmask++;
        }


        while (bitMask != 0) {
            int tIdx = idx + bitMask;
            bitMask >>= 1;
            if (tIdx > N)
                continue;
            if (cumFre >= tree[tIdx]) {
                // if the current cumulative frequency is equal to cumFre,
                // we are still looking for a higher index (if exists)
                idx = tIdx;
                cumFre -= tree[tIdx];
            }
        }
        if (cumFre != 0)
            return -1;
        else
            return idx;
    }


    // find the value of a particular array element
    int readSingle(int idx)
    {
        int sum = tree[idx]; // this sum will be decreased
        if (idx > 0) { // the special case
            int z = idx - (idx & -idx);
            idx--; // idx is not important anymore, so instead y, you can use idx
            while (idx != z) { // at some iteration idx (y) will become z
                sum -= tree[idx];
                // substruct tree frequency which is between y and "the same path"
                idx -= (idx & -idx);
            }
        }
        return sum;
    }

};


struct data {
    int x;
};

struct SEGTREE { // for range max with range addition update

    int N;
    vector<data> tree;
    vector<int> lazy;

    SEGTREE(int n)
    {
        N = n;
        tree.clear();
        tree.resize(4 * N + 3);
        lazy.clear();
        lazy.resize(4 * N + 3);
    }

    data make_data(int val)
    {
        data res;
        res.x = val;
        return res;
    }

    data combine(data a,data b)
    {
        data res;
        res.x = max(a.x, b.x);
        return res;
    }


    void push(int v) 
    {
        tree[v * 2].x += lazy[v];
        lazy[v * 2] += lazy[v];
        tree[v * 2 + 1].x += lazy[v];
        lazy[v * 2 + 1] += lazy[v];
        lazy[v] = 0;
    }

    void build(vi a, int v, int tl, int tr) 
    {
        if (tl == tr) {
            tree[v] = make_data(a[tl]);
        } else {
            int tm = (tl + tr) / 2;
            build(a, v * 2, tl, tm);
            build(a, v * 2 + 1, tm + 1, tr);
            tree[v] = combine(tree[v * 2], tree[v * 2 + 1]);
        }
    }


    void update(int v, int tl, int tr, int l, int r, int val) 
    {
        if (l > r)
            return;
        if (l == tl && tr == r) {
            tree[v].x += val;
            lazy[v] += val;
        } else {
            push(v);
            int tm = (tl + tr) / 2;
            update(v * 2, tl, tm, l, min(r, tm), val);
            update(v * 2 + 1, tm + 1, tr, max(l, tm + 1), r, val);
            tree[v] = combine(tree[v * 2], tree[v * 2 + 1]);
        }
    }

    // function overloading - same name, two update functions (one for range, and one for point)

    void update(int v, int tl, int tr, int pos, int new_val) 
    {
        if (tl == tr) {
            tree[v] = make_data(new_val);
        } else {
            int tm = (tl + tr) / 2;
            if (pos <= tm)
                update(v * 2, tl, tm, pos, new_val);
            else
                update(v * 2 + 1, tm + 1, tr, pos, new_val);
            tree[v] = combine(tree[v * 2], tree[v * 2 + 1]);
        }
    }

    data query(int v, int tl, int tr, int l, int r) 
    {
        if (l > r)
            return make_data(-inf);
        if (l <= tl && tr <= r)
            return tree[v];
        push(v);
        int tm = (tl + tr) / 2;
        return combine(query(v * 2, tl, tm, l, min(r, tm)),
                       query(v * 2 + 1, tm + 1, tr, max(l, tm + 1), r));
    }

};









//  simple square root decomposition for sum of subsegment

int n;
vector<int> a (n);

// preprocessing
int len = (int) sqrt (n + .0) + 1; // size of the block and the number of blocks
vector<int> b (len);
for (int i = 0; i < n; ++i)
    b[i / len] += a[i];

// answering the queries
for (;;) {
    int l, r; // input l and r
    int sum = 0;
    int c_l = l / len,   c_r = r / len;
    if (c_l == c_r)
        for (int i = l; i <= r; ++i)
            sum += a[i];
    else {
        for (int i = l, end = (c_l + 1) * len - 1; i <= end; ++i)
            sum += a[i];
        for (int i = c_l + 1; i <= c_r - 1; ++i)
            sum += b[i];
        for (int i = c_r * len; i <= r; ++i)
            sum += a[i];
    }
}

// Mo's algorithm simple implementation

void remove(idx);  // TODO: remove value at idx from data structure
void add(idx);     // TODO: add value at idx from data structure
int get_answer();  // TODO: extract the current answer of the data structure

int block_size;

struct Query {
    int l, r, idx;
    bool operator<(Query other) const
    {
        return make_pair(l / block_size, r) <
               make_pair(other.l / block_size, other.r);
    }
};

vector<int> mo_s_algorithm(vector<Query> queries) {
    vector<int> answers(queries.size());
    sort(queries.begin(), queries.end());

    // TODO: initialize data structure

    int cur_l = 0;
    int cur_r = -1;
    // invariant: data structure will always reflect the range [cur_l, cur_r]
    for (Query q : queries) {
        while (cur_l > q.l) {
            cur_l--;
            add(cur_l);
        }
        while (cur_r < q.r) {
            cur_r++;
            add(cur_r);
        }
        while (cur_l < q.l) {
            remove(cur_l);
            cur_l++;
        }
        while (cur_r > q.r) {
            remove(cur_r);
            cur_r--;
        }
        answers[q.idx] = get_answer();
    }
    return answers;
}

// faster sorting function for Mo's algorithm

bool cmp(pair<int, int> p, pair<int, int> q) {
    if (p.first / BLOCK_SIZE != q.first / BLOCK_SIZE)
        return p < q;
    return (p.first / BLOCK_SIZE & 1) ? (p.second < q.second) : (p.second > q.second);
}


struct SPARSE
{
    int N;
    int K;
    vector<int> lg;
    vector<vector<int>> st;

    int func(int val1, int val2)
    {
        int valr = max(val1, val2);
        return valr;
    }

    SPARSE(int n)
    {
        N = n;
        K = log2(N + 1);
        st.resize(N);
        for (int i = 0; i < N; i++)
        {
            st[i].resize(K + 1);
        }
        lg.resize(N + 1);
    }

    void build(vector<int> a)
    {
        for (int i = 0; i < N; i++)
            st[i][0] = a[i];

        for (int j = 1; j <= K; j++)
            for (int i = 0; i + (1 << j) <= N; i++)
                st[i][j] = func(st[i][j - 1], st[i + (1 << (j - 1))][j - 1]);

        lg[1] = 0;
        for (int i = 2; i <= N; i++)
            lg[i] = lg[i / 2] + 1;
    }

    int get(int l, int r)
    {
        int sum = 0;
        for (int j = K; j >= 0; j--) {
            if ((1 << j) <= r - l + 1) {
                sum = func(sum, st[l][j]);
                l += (1 << j);
            }
        }
    }


    // only works for idempotent functions
    int getFast(int l, int r)
    {
        int j = lg[r - l + 1];
        int ans = func(st[l][j], st[r - (1 << j) + 1][j]);
        return ans;
    }



};
// Disjoint set union using path compression


int parent[N]
int siz[N];
int rank[N];

int findLeader(int i)
{
    if (parent[parent[i]] != parent[i])
        parent[i] = findLeader(parent[i]);
    return parent[i];
}

// union by size

void unionSets(int a, int b)
{
    int parent_a = findLeader(a), parent_b = findLeader(b);
    if (parent_a == parent_b)
        return;
    if (siz[parent_a] >= siz[parent_b]) {
        swap(parent_a, parent_b);
    }
    siz[parent_b] += siz[parent_a];
    parent[parent_a] = parent_b;
    return;
}

// union by rank

void unionSets(int a, int b)
{
    int parent_a = findLeader(a), parent_b = findLeader(b);
    if (parent_a == parent_b)
        return;
    if (rank[parent_a] >= rank[parent_b]) {
        swap(parent_a, parent_b);
    }

    if (rank[parent_b] == rank[parent_a])
        rank[parent_b]++;
    parent[parent_a] = parent_b;
    return;
}

// kruskals algorithm using dsu for minimum spanning tree

vector<int> parent, rank;

void make_set(int v) {
    parent[v] = v;
    rank[v] = 0;
}

int find_set(int v) {
    if (v == parent[v])
        return v;
    return parent[v] = find_set(parent[v]);
}

void union_sets(int a, int b) {
    a = find_set(a);
    b = find_set(b);
    if (a != b) {
        if (rank[a] < rank[b])
            swap(a, b);
        parent[b] = a;
        if (rank[a] == rank[b])
            rank[a]++;
    }
}

struct Edge {
    int u, v, weight;
    bool operator<(Edge const& other) {
        return weight < other.weight;
    }
};

int n;
vector<Edge> edges;

int cost = 0;
vector<Edge> result;
parent.resize(n);
rank.resize(n);
for (int i = 0; i < n; i++)
    make_set(i);

sort(edges.begin(), edges.end());

for (Edge e : edges) {
    if (find_set(e.u) != find_set(e.v)) {
        cost += e.weight;
        result.push_back(e);
        union_sets(e.u, e.v);
    }
}

// Dijkstra using Priority Queue
typedef pair<int, int> iPair;

// To add an edge
void addEdge(vector <pair<int, int> > adj[], int u,
             int v, int wt)
{
    adj[u].push_back(make_pair(v, wt));
    adj[v].push_back(make_pair(u, wt));
}


// Prints shortest paths from src to all other vertices
void shortestPath(vector<pair<int, int> > adj[], int V, int src)
{
    // Create a priority queue to store vertices that
    // are being preprocessed. This is weird syntax in C++.
    // Refer below link for details of this syntax
    // http://geeksquiz.com/implement-min-heap-using-stl/
    priority_queue< iPair, vector <iPair> , greater<iPair> > pq;

    // Create a vector for distances and initialize all
    // distances as infinite (INF)
    vector<int> dist(V, INF);

    // Insert source itself in priority queue and initialize
    // its distance as 0.
    pq.push(make_pair(0, src));
    dist[src] = 0;

    /* Looping till priority queue becomes empty (or all
    distances are not finalized) */
    while (!pq.empty())
    {
        // The first vertex in pair is the minimum distance
        // vertex, extract it from priority queue.
        // vertex label is stored in second of pair (it
        // has to be done this way to keep the vertices
        // sorted distance (distance must be first item
        // in pair)
        int u = pq.top().second;
        pq.pop();

        // Get all adjacent of u.
        for (auto x : adj[u])
        {
            // Get vertex label and weight of current adjacent
            // of u.
            int v = x.first;
            int weight = x.second;

            // If there is shorted path to v through u.
            if (dist[v] > dist[u] + weight)
            {
                // Updating distance of v
                dist[v] = dist[u] + weight;
                pq.push(make_pair(dist[v], v));
            }
        }
    }

    // Print shortest distances stored in dist[]
    printf("Vertex Distance from Source\n");
    for (int i = 0; i < V; ++i)
        printf("%d \t\t %d\n", i, dist[i]);
}


// Breadth first search

vector<vector<int>> adj;  // adjacency list representation
int n; // number of nodes
int s; // source vertex

queue<int> q;
vector<bool> used(n);
vector<int> d(n), p(n);

q.push(s);
used[s] = true;
p[s] = -1;
while (!q.empty()) {
    int v = q.front();
    q.pop();
    for (int u : adj[v]) {
        if (!used[u]) {
            used[u] = true;
            q.push(u);
            d[u] = d[v] + 1;
            p[u] = v;
        }
    }
}

// path printing in bfs

if (!used[u]) {
    cout << "No path!";
} else {
    vector<int> path;
    for (int v = u; v != -1; v = p[v])
        path.push_back(v);
    reverse(path.begin(), path.end());
    cout << "Path: ";
    for (int v : path)
        cout << v << " ";
}



//  Deapth first search (generic implementation)

vector<vector<int>> adj; // graph represented as an adjacency list
int n; // number of vertices

vector<int> color;

vector<int> time_in, time_out;
int dfs_timer = 0;

void dfs(int v) {
    time_in[v] = dfs_timer++;
    color[v] = 1;
    for (int u : adj[v])
        if (color[u] == 0)
            dfs(u);
    color[v] = 2;
    time_out[v] = dfs_timer++;
}

