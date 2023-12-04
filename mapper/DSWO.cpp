#include <stdio.h>
#include <vector>
#include <queue>
#include <stack>
#include <map>
#include <set>
#include <algorithm>
#include <string.h>
#include <iostream>
#include <sstream>
#include <cassert>
#include <random>
#include <iterator>
#include <string>
#include <ctime>
#include <chrono>
using namespace std;

const double sample_eps = 1e-10;

const int MAX_QUBIT = 1225;
const int MAX_GATES = 1440000;

//decay coeff
const double decayCoeff = 0.001;
const int decayCycle = 5;

//Dist
double D[MAX_QUBIT][MAX_QUBIT];

//LINK STATE
bool LAYOUT[MAX_QUBIT][MAX_QUBIT];

//EXTEND_SIZE
int SIZEE = 20;

//DSWO parameter
int LIMIT_COUNT = 15;


/*
    Swap struct
*/
struct Swap{
    int u, v;
    Swap(int _u, int _v){
        u = _u;
        v = _v;
    }
    Swap(){u = 0;v = 0;}
};

/*
    Swap with score
*/
struct TempSwap {
    double heuScore;
    Swap swap;

    TempSwap(){}
};

/*
    EDGES
*/
vector<Swap> BIGEDGES[MAX_QUBIT];

/*
    GATE
*/
struct Gate{
    /* data */
    int u, v, id;
    int qubitNumber;
    short preNumber;
    vector<int> edges;
    vector<int> attach;
    Gate(int _u){
        u = _u;
        qubitNumber = 1;
    }
    Gate(int _u, int _v){
        u = _u;
        v = _v;
        qubitNumber = 2;
    }
    Gate(){}
};


map<int, Gate> GATES;

/*
    Gate
*/
struct AnsGate{
    int u, v, id;
    int qubitNumber;
    AnsGate(int _u){
        u = _u;
        v = -1;
        qubitNumber = 1;
    }
    AnsGate(int _u, int _v){
        u = _u;
        v = _v;
        if (_v == -1){
            qubitNumber = 1;
        }else{
            qubitNumber = 2;
        }
    }

    AnsGate(const AnsGate &target){
        u = target.u;
        v = target.v;
        qubitNumber = target.qubitNumber;
        id = target.id;
    }
};

vector<AnsGate> initialGates;
vector<AnsGate> ansTrees;

inline void addAnsGate(AnsGate &ansG){
    ansTrees.push_back(ansG);
}

/*
    dealH优化
*/
struct WeightGate{
    int gateid;
    double weight;
    WeightGate(int _gateid, double _weight){
        gateid = _gateid;
        weight = _weight;
    }
};

/*
    SABRE STATE
*/

struct SABRESTATE{
    /*
        Front Layer的id
    */
    vector<int> Fids;
    int ansf1;

    //initial gates, qubit number
    int preCount, qubitNumber;

    //topology mapping
    int l2p[MAX_QUBIT], p2l[MAX_QUBIT];

    int breakTime;
    int breakAllCount;
    int breakAllGCount;

    //inDegree
    short preNumber[MAX_GATES];

    //execution gates, deadloop count
    int count2, countSwap, countDeadLoop;

    //heuScore
    double heuScore;

    //Extend Set
    vector<int> ESet;
    
    //optGate
    vector<WeightGate> readyGates[MAX_QUBIT];

    //preDealH
    double preDealH;

    //decay var
    double decay[MAX_QUBIT];
    int depth[MAX_QUBIT];
    int nowDepth;
    int decayTime;

    void printl2p(){
        for (int i = 0;i < qubitNumber;++i){
            printf("%d ", l2p[i]);
        }
        printf("\n");
        for (int i = 0;i < qubitNumber;++i){
            printf("%d ", p2l[i]);
        }
        printf("\n");
    }

    SABRESTATE(int _qubitNumber, vector<int> _Fids, int * _l2p, short *_preNumber, int _preCount){
        qubitNumber = _qubitNumber;
        preCount = _preCount;

        ansf1 = 0;
        for (int i = 0;i < qubitNumber;++i){
            decay[i] = 1;
            depth[i] = 0;
        }
        nowDepth = 0;

        for (int i = 0;i < qubitNumber;++i){
            l2p[i] = _l2p[i];
        }
        for (int i = 0;i < qubitNumber;++i){
            p2l[l2p[i]] = i;
        }

        memcpy(preNumber, _preNumber, preCount * sizeof(short));
        Fids = _Fids;

        count2 = 0;
        countSwap = 0;
        decayTime = 0;
        countDeadLoop = 0;
        breakAllCount = 0;
        breakAllGCount = 0;
        breakTime = 0;
    }

    bool canExecute(Gate gate){
        return LAYOUT[l2p[gate.u]][l2p[gate.v]];
    }

    void preDeal(int sizeE){

        for (int i = 0;i < qubitNumber;++i){
            readyGates[i].clear();
        }

        double preH = 0;
        int FCount = Fids.size();
        queue<int> EQueue;
        for (auto gateid:Fids){
            Gate gate = GATES[gateid];
            preH += D[l2p[gate.u]][l2p[gate.v]] / FCount;
            WeightGate wg = WeightGate(gateid, 1.0 / FCount);
            readyGates[l2p[gate.u]].push_back(wg);
            readyGates[l2p[gate.v]].push_back(wg);
            EQueue.push(gateid);
        }

        ESet.clear();
        vector<int> decQueue;

        while (ESet.size() < sizeE && EQueue.size()){
            int gateid = EQueue.front();
            Gate gate = GATES[gateid];
            EQueue.pop();
            decQueue.push_back(gateid);
            for (auto succid:gate.edges){
                Gate succ = GATES[succid];
                --preNumber[succ.id];
                if (preNumber[succ.id] < 0){
                    printf("--%d\n", succid);
                }
                assert(preNumber[succ.id] >= 0);
                if (preNumber[succ.id] == 0){
                    ESet.push_back(succ.id);
                    EQueue.push(succ.id);
                }
            }
        }

        if (ESet.size() > sizeE){
            ESet = vector<int>({ESet.begin(), ESet.begin() + sizeE});
        }

        int ECount = ESet.size();
        for (auto gateid:ESet){
            Gate gate = GATES[gateid];
            preH += D[l2p[gate.u]][l2p[gate.v]] / ECount;
            WeightGate wg = WeightGate(gateid, 0.5 / ECount);
            //printf("%d %d %f %f\n", l2p[gate.u], l2p[gate.v], wg.weight, preH);
            readyGates[l2p[gate.u]].push_back(wg);
            readyGates[l2p[gate.v]].push_back(wg);
        }

        for (auto gateid:decQueue){
            Gate gate = GATES[gateid];
            for (auto gateidd:gate.edges){
                //printf("++%d \n", gateidd);
                ++preNumber[gateidd];
                //printf("delete id %d", gateidd)
            }
        }
        preDealH = preH;
        //printf("preDealH: %f\n", preDealH);
    }

    double preWithH(Swap swap){
        double dealH = preDealH;
        int u = swap.u;
        int v = swap.v;
        //printf("%d %d %f\n", u, v, preDealH);x
        for (auto &wg:readyGates[u]){
            Gate gate = GATES[wg.gateid];
            double coff = wg.weight;
            int vv = l2p[gate.u];
            if (vv == u)
                vv = l2p[gate.v];
            if (vv == v)
                continue;
            // [u][vv] => [v][vv]

            //printf("ready: %d %d %f\n", gate.u, gate.v, wg.weight);
            dealH += coff * (D[v][vv] - D[u][vv]);
        }
        for (auto &wg:readyGates[v]){
            Gate gate = GATES[wg.gateid];

            double coff = wg.weight;
            int uu = l2p[gate.u];
            if (uu == v)
                uu = l2p[gate.v];
            if (uu == u)
                continue;
            // [v][uu] => [u][uu]

            //printf("ready: %d %d %f\n", gate.u, gate.v, wg.weight);
            dealH += coff * (D[u][uu] - D[v][uu]);
        }

        return dealH;
    }

    vector<Swap> obtainSwaps(){
        vector<Swap> candidates;
        set<int> bits;
        for (auto gateid:Fids){
            Gate gate = GATES[gateid];
            bits.insert(l2p[gate.u]);
            bits.insert(l2p[gate.v]);
        }
        for (auto id:bits){
            for (auto swap:BIGEDGES[id]){
                int v = swap.v;
                int u = swap.u;
                if (bits.find(v) != bits.end() && bits.find(u) != bits.end() && u > v){

                    continue;
                }
                candidates.push_back(swap);
            }
        }
        return candidates;
    }

    void executeSwap(Swap theSwap){
        int u = theSwap.u;
        int v = theSwap.v;

        l2p[p2l[u]] = v;
        l2p[p2l[v]] = u;

        int t = p2l[u];
        p2l[u] = p2l[v];
        p2l[v] = t;

        AnsGate ansGate = AnsGate(u, v);
        ansGate.id = -1;

        addAnsGate(ansGate);
        ansf1 += 1;

        countSwap += 1;
        //printl2p();
    }

    set<int> cancelGate(queue<AnsGate> &gates){
        set<int> allowG;
        allowG.clear();

        //int preL2p[MAX_QUBIT];
        //memcpy(preL2p, l2p, qubitNumber * sizeof(int));

        int size = gates.size();
        /*if (size == 5){
            printf("cancel ??F:\n");
            for (auto fid:Fids){
                Gate gate = GATES[fid];
                printf("%d %d | ", l2p[gate.u], l2p[gate.v]);
            }
        }*/

        while (!gates.empty()){
            auto gate = gates.front();
            allowG.insert(gate.u);
            allowG.insert(gate.v);

            gates.pop();

            //printl2p();

            int u = gate.u;
            int v = gate.v;

            /*if (size == 5){
                printf("cancel: %d %d\n", u, v);
            }*/

            l2p[p2l[u]] = v;
            l2p[p2l[v]] = u;

            int t = p2l[u];
            p2l[u] = p2l[v];
            p2l[v] = t;

            //printl2p();
        }
        /*for(int i = 0;i < qubitNumber;++i){
            out += D[l2p[i]][preL2p[i]];
        }*/
        return allowG;
    }

    bool findShortestPath(int size, set<int> &allowG){
        ansf1 = 0;
        int shortest = qubitNumber * qubitNumber;
        int u = 0, v = 0;
        for (auto fid:Fids){
            auto gate = GATES[fid];
            /*if (allowG.count(l2p[gate.u]) == 0 && allowG.count(l2p[gate.v]) == 0){
                continue;
            }*/

            if (D[l2p[gate.u]][l2p[gate.v]] < shortest){
                u = l2p[gate.u];
                v = l2p[gate.v];
                shortest = D[l2p[gate.u]][l2p[gate.v]];
            }
        }
        assert(u != v);

        int DD[MAX_QUBIT + 1];
        int PRE[MAX_QUBIT + 1];
        for (int i = 0;i < qubitNumber;++i){
            DD[i] = qubitNumber * qubitNumber;
            PRE[i] = -1;
        }

        queue<int> queues;
        while (!queues.empty()){
            queues.pop();
        }
        queues.push(u);
        DD[u] = 0;
        while (queues.size() > 0){
            int head = queues.front();
            queues.pop();

            for (auto &swap:BIGEDGES[head]){
                auto goal = swap.v;
                if (PRE[goal] == -1){
                    queues.push(goal);
                }
                if (DD[head] + 1 < DD[goal]){
                    DD[goal] = DD[head] + 1;
                    PRE[goal] = head;
                }

                if (goal == v){
                    break;
                }
            }
            if (PRE[v] != -1){
                break;
            }
        }
        auto now = v;

        vector<int> path;
        path.clear();
        path.push_back(now);
        while (PRE[now] != u){
            auto goal = PRE[now];
            path.push_back(goal);
            now = goal;
        }
        path.push_back(u);
        int i = 0;
        int j = path.size();

        /*if (size == 5){
            printf("new path:%d\n", j);
            for (int _i = 0;_i < j;++_i){
                printf("path:%d\n", path[_i]);
            }
        }*/
        --j;
        while (i + 1 < j){
            auto swap = Swap(path[i], path[i+1]);
            executeSwap(swap);
            ++i;
            if (i + 1 >= j){
                break;
            }
            swap = Swap(path[j], path[j-1]);
            executeSwap(swap);
            --j;
        }

        auto nowCount = count2;
        executeGates();
        auto newCount = count2;
        assert(newCount > nowCount);
        return true;
    }

    void checkDeadLoop(){
        if (ansf1 >= 10 * qubitNumber){
            int lastIndex = ansTrees.size() - 1;
            auto nowGate = ansTrees[lastIndex--];
            if (nowGate.id != -1){
                return;
            }

            queue<AnsGate> cancelGates;
            while (!cancelGates.empty()){
                cancelGates.pop();
            };
            cancelGates.push(nowGate);

            while (lastIndex >= 0){
                auto newGateNode = ansTrees[lastIndex];

                if (newGateNode.id != -1){
                    break;
                }
                --lastIndex;
                cancelGates.push(newGateNode);
            }
            ansTrees.erase(ansTrees.begin() + lastIndex + 1, ansTrees.end());
            countDeadLoop += 1;
            breakTime += 1;

            int size = cancelGates.size();

            set<int> allowG = cancelGate(cancelGates);

            breakAllGCount += size;

            findShortestPath(size, allowG);

            //printf("deadloop end\n");

            return;
        }

        if (ansTrees.size() == 0){
            return;
        }
        int lastIndex = ansTrees.size() - 1;
        auto nowGate = ansTrees[lastIndex--];
        //printf("try deadloopcheck1 %d %d %d %ld\n", nowGate.id, nowGate.u, nowGate.v, ansTrees.size());
        if (nowGate.id != -1){
            return;
        }

        queue<AnsGate> cancelGates;
        while (!cancelGates.empty()){
            cancelGates.pop();
        };
        cancelGates.push(nowGate);

        int count = abs(LIMIT_COUNT);

        while (lastIndex >= 0 && count > 0){
            auto newGateNode = ansTrees[lastIndex];
            if (newGateNode.id != -1){
                return;
            }
            cancelGates.push(newGateNode);
            if ((newGateNode.u == nowGate.u && newGateNode.v == nowGate.v) || (newGateNode.u == nowGate.v && newGateNode.v == nowGate.u)){
                countDeadLoop += 1;
                ansTrees.erase(ansTrees.begin() + lastIndex, ansTrees.end());
                breakTime += 1;

                int size = cancelGates.size();

                breakAllGCount += size;

                set<int> allowG = cancelGate(cancelGates);

                //breakAllCount += newCount;
                //printf("break count:%d %d\n", newCount);
                findShortestPath(size, allowG);
                return;
            }
            count -= 1;
            lastIndex -= 1;
        }
    }

    TempSwap generateNextTempSwap(Swap swap){
        TempSwap tempSwap = TempSwap();
        tempSwap.swap = swap;
        heuScore = preWithH(swap) * max(decay[swap.u], decay[swap.v]);
        if (heuScore < 0){
            assert(heuScore >= 0);
        }

        tempSwap.heuScore = heuScore;
        return tempSwap;
    }

    void executeTempSwap(Swap swap){
        executeSwap(swap);
        decay[swap.u] += decayCoeff;
        decay[swap.v] += decayCoeff;
        decayTime += 1;

        int nowDepth = max(depth[p2l[swap.u]], depth[p2l[swap.v]]) + 1;
        depth[p2l[swap.u]] = depth[p2l[swap.v]] = nowDepth;
        nowDepth = max(nowDepth, nowDepth);

        if (decayTime % decayCycle == 0){
            for (int i = 0;i < qubitNumber;++i){
                decay[i] = 1;
            }
        }
        checkDeadLoop();
    }

    vector<TempSwap> searchNextSwaps(){
        vector<TempSwap> temps;

        vector<Swap> swaps = obtainSwaps();
        preDeal(SIZEE);
        for (auto &swap:swaps){
            auto temp = generateNextTempSwap(swap);
            temps.push_back(temp);
        }

        return temps;
    } 

    void executeGates(){

        vector<int> exeGateList;
        while (1){
            exeGateList.clear();

            for (auto gateid:Fids){
                Gate gate = GATES[gateid];

                if (canExecute(gate)){
                    exeGateList.push_back(gateid);
                    AnsGate ansGate = AnsGate(l2p[gate.u], l2p[gate.v]);
                    ansGate.id = gate.id;

                    addAnsGate(ansGate);

                    count2 += 1;
                    for (auto sgid:gate.attach){
                        Gate sg = GATES[sgid];
                        AnsGate ansGate = AnsGate(l2p[sg.u]);
                        ansGate.id = sg.id;
                        addAnsGate(ansGate);
                    }
                }
            }

            if (exeGateList.size() == 0)
                break;

            ansf1 = 0;
            for (auto gid:exeGateList){

                Fids.erase(remove(Fids.begin(), Fids.end(), gid), Fids.end());
                Gate g = GATES[gid];
                for (auto succid:g.edges){
                    Gate succ = GATES[succid];
                    --preNumber[succ.id];
                    assert(preNumber[succ.id] >= 0);
                    if (preNumber[succ.id] == 0){
                        Fids.push_back(succ.id);
                    }
                }
            }
            for (int i = 0;i < qubitNumber;++i){
                decay[i] = 1;
            }
        }
        /*while (Fids.size() > 0 && findShortestPath(0)){

        }*/
    }
};

bool comp(const TempSwap &a, const TempSwap &b){
    return a.heuScore < b.heuScore;
}

void searchState(SABRESTATE &state){
    if (state.Fids.size() == 0){
        return;
    }
    int count = 0;
    while (1){
        count += 1;
        vector<TempSwap> queues = state.searchNextSwaps();
        sort(queues.begin(), queues.end(), comp);

        assert(queues.size() > 0);
        double score = queues[0].heuScore;
        int i = 0;
        for (auto &tSwap:queues){
            double newScore = tSwap.heuScore;
            if (abs(newScore - score) > sample_eps)
                break;
            i += 1;
        }
        //cout<<queues.size()<<endl;

        TempSwap tswap;

        if (i > 1){
            tswap = queues[rand() % i];
        }else{
            tswap = queues[0];
        }

        state.executeTempSwap(tswap.swap);

        state.executeGates();
        if (state.Fids.size() == 0){
            return;
        }
    }
}

extern "C"
{
    void freeDSWO(char* ptr){
        free(ptr);
    }
}

extern "C"
{   
    char * DSWO(int *gateContent, bool *layout, int *l2p, int gateLen, int qubitNumber, int limitL){
        //printf("in acce SABRE %d %d %d\n", gateLen, qubitNumber, limitL);

        srand((unsigned)time(0));
        setbuf(stdout, NULL);
        LIMIT_COUNT = limitL;

        ansTrees.clear();

        for (int i = 0;i < qubitNumber;++i){
            BIGEDGES[i].clear();
            for (int j = 0;j < qubitNumber;++j){
                LAYOUT[i][j] = layout[i * qubitNumber + j];
                D[i][j] = LAYOUT[i][j] ? 1 : qubitNumber * qubitNumber * 2;
                if (i == j){
                    D[i][j] = 0;
                }
                if (LAYOUT[i][j] && i != j){
                    Swap swap = Swap(i, j);
                    BIGEDGES[i].push_back(swap);
                }
            }
        }

        for (int k = 0;k < qubitNumber;++k){
            for (int i = 0;i < qubitNumber;++i){
                if (i == k)
                    continue;
                for (int j = 0;j < qubitNumber;++j){
                    if (j == k || j == i){
                        continue;
                    }
                    D[i][j] = min(D[i][j], D[i][k] + D[k][j]);
                }
            }
        }

        vector<AnsGate> ansGates;
        ansGates.clear();

        for (int i = 0;i < gateLen;++i){
            int u, v, id;
            u = gateContent[3 * i];
            v = gateContent[3 * i + 1];
            id = gateContent[3 * i + 2];
            AnsGate ansGate = AnsGate(u, v);
            ansGate.id = id;
            ansGates.push_back(ansGate);
        }

        int predag[MAX_QUBIT];
        initialGates.clear();
        for (int i = 0;i < qubitNumber;++i){
            predag[i] = 0;
        }

        short preNumbers[MAX_GATES];

        vector<int> Fids;
        Fids.clear();
        GATES.clear();

        for (auto &aG:ansGates){
            Gate gg = Gate();
            gg.u = aG.u;
            gg.v = aG.v;
            gg.id = aG.id;
            gg.qubitNumber = aG.qubitNumber;
            gg.edges = vector<int>();
            gg.attach = vector<int>();

            GATES[gg.id] = gg;
            short preNumber = 0;
            if (gg.qubitNumber == 1){
                int gateid = predag[gg.u];
                if (gateid != 0){
                    (&GATES[gateid])->attach.push_back(gg.id);
                }else{
                    aG.u = l2p[aG.u];
                    initialGates.push_back(aG);
                }
            }else{
                int qubits[2] = {aG.u, aG.v};
                for (int i = 0;i < 2;++i){
                    int qubit = qubits[i];
                    int gateid = predag[qubit];
                    if (gateid != 0){
                        auto pos = find(GATES[gateid].edges.begin(), GATES[gateid].edges.end(), gg.id);
                        if (pos == (&GATES[gateid])->edges.end()){
                            //
                            (&GATES[gateid])->edges.push_back(gg.id);
                            preNumber += 1;
                        }
                    }
                }

                predag[aG.u] = gg.id;
                predag[aG.v] = gg.id;
                gg.preNumber = preNumber;

                if (preNumber == 0){
                    Fids.push_back(gg.id);
                }
            }
            preNumbers[gg.id] = preNumber;
        }
        SABRESTATE initialSABRE = SABRESTATE(qubitNumber, Fids, l2p, preNumbers, gateLen + 1);
        initialSABRE.executeGates();
        searchState(initialSABRE);

        //printf("average count:%d %d %d\n", initialSABRE.breakAllGCount, initialSABRE.breakAllCount, initialSABRE.breakTime);

        for (int i = 0;i < ansTrees.size();++i){
            initialGates.push_back(AnsGate(ansTrees[i]));
        }

        stringstream fmt;
        fmt.str("");
        fmt<<initialGates.size()<<".";
        for (auto &ag:initialGates){
            if (ag.id == -1){
                fmt<<ag.u<<","<<ag.v<<".";
            }else{
                fmt<<ag.id<<","<<ag.u<<","<<ag.v<<".";
            }
        }
        string output;
        fmt>>output;
        //cout<<output<<endl;
        int len = output.length();
        char *store = (char *)malloc((len+1)*sizeof(char));
        //output.copy(store, len, 0);
        //char *new_buf = strdup(store);
        //printf("allocated address: %p\n", new_buf);
        strcpy(store, output.c_str());

        return store;
    }

}
