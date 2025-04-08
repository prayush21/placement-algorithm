/***************************************************************
 * Fiduccia–Mattheyses (FM) Partitioning with:
 *  - Two Area Modes (1 => each gate area=1, 2 => width*height)
 *  - Optional ML Heuristic (logistic gating).
 *  - Multi-Pass Iteration.
 *
 * * Compile Example:
 *   g++ -std=c++17 partitioner.cpp -o partitioner -O2
 *
 * Usage:
 *   ./partitioner <benchmarkName> <areaMode=1|2> <mlMode=0|1> <maxPasses>
 *
 * Example:
 *   ./partitioner superblue18 2 1 5
 *

 * Code Flow:
 *   1) Reads <benchmark>/<benchmark>.nodes for cell info.
 *   2) Depending on "areaMode", sets each cell's area to 1 or (width*height).
 *   3) Creates an initial partition (random) with ~10% area imbalance allowed.
 *   4) Reads <benchmark>/<benchmark>.nets to build net->cells adjacency.
 *   5) Repeatedly runs an FM pass up to <maxPasses> times:
 *      - Each pass picks highest-gain cells (unlocked),
 *      - Checks area feasibility,
 *      - Optionally gates with an ML logistic model,
 *      - Updates neighbor gains.
 *      - Tracks the best cut and reverts if no improvement.
 *   6) Prints the final cut size and runtime.
 *
 ***************************************************************/

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <map>
#include <random>
#include <algorithm>
#include <chrono>
#include <limits>
#include <sstream>
#include <cmath>
#include <ctime>
#include <list>

// ----------------------------------------------------------
// Cell and Net Classes
// ----------------------------------------------------------
class Cell
{
public:
    std::string id;

    // area=1 if areaMode=1, or width*height if areaMode=2
    long long area;

    // FM data
    int cell_gain;   // computed gain
    int partition;   // 0 => Partition A, -1 => Partition B
    int lock_status; // 0 => unlocked, 1 => locked

    // For adjacency
    std::unordered_set<std::string> net_list;

    // Iterator to the gain bucket
    std::list<std::string>::iterator bucket_iterator;

    Cell()
        : area(1), cell_gain(0), partition(0), lock_status(0)
    {
    }

    Cell(const std::string &cid)
        : id(cid), area(1), cell_gain(0), partition(0), lock_status(0)
    {
    }
};

class Net
{
public:
    std::string id;

    // All cell IDs in this net
    std::vector<std::string> cell_list;

    // #cells in partition A or B
    int Asize;
    int Bsize;

    Net() : Asize(0), Bsize(0) {}
    Net(const std::string &netId)
        : id(netId), Asize(0), Bsize(0)
    {
    }
};

// ----------------------------------------------------------
// Global Data
// ----------------------------------------------------------

// Master maps
std::unordered_map<std::string, Cell> cell_map;
std::unordered_map<std::string, Net> net_map;

// Keep track of total area in each partition
long long areaA = 0;
long long areaB = 0;

// Balancing ratio (10% by default)
double ALLOWED_AREA_IMBALANCE_RATIO = 0.2;
int P_MAX = -1;
int max_gain_index = 0;
int min_gain_index = 0;

//  std::map<int, std::list<std::string>> gain_bucket;
std::vector<std::list<std::string>> gain_bucket;

// For multi-pass tracking
std::unordered_map<std::string, int> bestPartitionRecord; // cell->partition
int bestCutGlobal = std::numeric_limits<int>::max();
long long bestAreaA = 0;
long long bestAreaB = 0;

// Prototypes
void readCells(const std::string &nodefile, int areaMode);
void readNets(const std::string &netfile);
void initPartition();
void printPartition(const std::string &partitionName, bool showLockStatus);
void printGainBucket();
void printGainBucketMaxMin();
int computeCutSize();
bool canMoveCellBalanced(const std::string &cellId);
void initializeGainBuckets();
int runOnePass(int mlMode);
double mlScoreCell(const Cell &c);
void updateNeighborGains(const std::string &cellId, int oldPartition);
void updateGainsAfterMove(const std::string &base_cell_id);
void updateCellGain(Cell &cell, int delta);

// ----------------------------------------------------------
// 1) Reading .nodes => populates cell_map
// ----------------------------------------------------------
void readCells(const std::string &nodefile, int areaMode)
{
    cell_map.clear();
    areaA = 0;
    areaB = 0;

    std::ifstream fin(nodefile);
    if (!fin.is_open())
    {
        std::cerr << "[ERROR] Cannot open node file: " << nodefile << std::endl;
        exit(1);
    }

    std::string line;
    while (std::getline(fin, line))
    {
        if (line.empty() || line[0] == '#' ||
            line.find("Num") != std::string::npos ||
            line.find("UCLA") != std::string::npos)
        {
            continue;
        }

        // Format e.g. "cell0 width height"
        std::stringstream ss(line);
        std::string cellId;
        long long width, height;
        ss >> cellId >> width >> height;
        if (cellId.empty())
        {
            continue;
        }

        Cell c(cellId);
        // areaMode=1 => all area=1
        // areaMode=2 => area=width*height
        if (areaMode == 2)
        {
            long long computed = width * height;
            if (computed <= 0)
                computed = 1;
            c.area = computed;
        }
        else
        {
            c.area = 1; // fallback
        }

        c.cell_gain = 0;
        c.partition = 0;
        c.lock_status = 0;
        cell_map[cellId] = c;
    }
    fin.close();

    std::cout << "[INFO] Read " << cell_map.size()
              << " cells from " << nodefile
              << " with areaMode=" << areaMode << std::endl;
}

// ----------------------------------------------------------
// 2) Reading .nets => populates net_map & updates net->cell list
// ----------------------------------------------------------
void readNets(const std::string &netfile)
{
    net_map.clear();

    std::ifstream fin(netfile);
    if (!fin.is_open())
    {
        std::cerr << "[ERROR] Cannot open net file: " << netfile << std::endl;
        exit(1);
    }

    std::string line;
    std::string currentNetId;
    while (std::getline(fin, line))
    {
        if (line.empty() || line[0] == '#' ||
            line.find("Num") != std::string::npos ||
            line.find("UCLA") != std::string::npos)
        {
            continue;
        }

        std::stringstream ss(line);
        std::string token;
        ss >> token;

        if (token == "NetDegree")
        {
            // Format: "NetDegree : X net1234"
            std::string colon;
            int netDegree;
            ss >> colon >> netDegree >> currentNetId;
            Net netObj(currentNetId);
            netObj.Asize = 0;
            netObj.Bsize = 0;
            net_map[currentNetId] = netObj;
        }
        else
        {
            std::string cellId = token;
            if (cell_map.find(cellId) == cell_map.end())
            {
                // Possibly a fixed cell/pad not in cell_map
                continue;
            }
            // Add to net->cells
            net_map[currentNetId].cell_list.push_back(cellId);
            // Add net to cell->nets
            cell_map[cellId].net_list.insert(currentNetId);

            // Update the net's Asize/Bsize depending on partition
            if (cell_map[cellId].partition == 0)
            {
                net_map[currentNetId].Asize++;
            }
            else
            {
                net_map[currentNetId].Bsize++;
            }
        }
    }
    fin.close();

    std::cout << "[INFO] Read " << net_map.size()
              << " nets from " << netfile << std::endl;
}

// ----------------------------------------------------------
// 3) Initial Partition (random) while respecting ~(ALLOWED_AREA_IMBALANCE_RATIO*100)% imbalance
// ----------------------------------------------------------
void initPartition()
{
    areaA = 0;
    areaB = 0;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, 1);

    for (auto &cpair : cell_map)
    {
        auto &c = cpair.second;
        int side = dis(gen); // 0 or 1

        if (side == 0)
        {
            long long newA = areaA + c.area;
            long long diff = std::llabs(newA - areaB);
            long long total = newA + areaB;
            if (total == 0)
            {
                // trivial
                c.partition = 0;
                areaA = newA;
            }
            else
            {
                double ratio = (double)diff / (double)total;
                if (ratio <= ALLOWED_AREA_IMBALANCE_RATIO)
                {
                    c.partition = 0;
                    areaA = newA;
                }
                else
                {
                    c.partition = -1;
                    areaB += c.area;
                }
            }
        }
        else
        {
            long long newB = areaB + c.area;
            long long diff = std::llabs(areaA - newB);
            long long total = areaA + newB;
            if (total == 0)
            {
                c.partition = -1;
                areaB = newB;
            }
            else
            {
                double ratio = (double)diff / (double)total;
                if (ratio <= ALLOWED_AREA_IMBALANCE_RATIO)
                {
                    c.partition = -1;
                    areaB = newB;
                }
                else
                {
                    c.partition = 0;
                    areaA += c.area;
                }
            }
        }
        c.lock_status = 0;
    }
}

// ----------------------------------------------------------
// 4) Compute cut size: a net is cut if net.Asize>0 & net.Bsize>0
// ----------------------------------------------------------
int computeCutSize()
{
    int cut = 0;
    for (const auto &n : net_map)
    {
        if (n.second.Asize > 0 && n.second.Bsize > 0)
        {
            cut++;
        }
    }
    return cut;
}

// ----------------------------------------------------------
// 5) Check if we can move a cell while respecting ~(ALLOWED_AREA_IMBALANCE_RATIO*100)% imbalance
//    If feasible, we update areaA/areaB
// ----------------------------------------------------------
bool canMoveCellBalanced(const std::string &cellId)
{
    auto &c = cell_map[cellId];
    long long newA = areaA;
    long long newB = areaB;

    if (c.partition == 0)
    {
        // Move from A->B
        newA -= c.area;
        newB += c.area;
    }
    else
    {
        // Move from B->A
        newB -= c.area;
        newA += c.area;
    }

    if (newA < 0 || newB < 0)
    {
        return false;
    }

    long long total = newA + newB;
    if (total == 0)
    {
        // trivial
        areaA = newA;
        areaB = newB;
        return true;
    }
    long long diff = std::llabs(newA - newB);
    double ratio = (double)diff / (double)total;

    if (ratio <= ALLOWED_AREA_IMBALANCE_RATIO)
    {
        // update global
        areaA = newA;
        areaB = newB;
        return true;
    }
    return false;
}

// ----------------------------------------------------------
// 6) Initialize Gains + Build Buckets
// ----------------------------------------------------------
void initializeGainBuckets()
{
    gain_bucket.clear();

    // Re-init net Asize/Bsize
    for (auto &np : net_map)
    {
        np.second.Asize = 0;
        np.second.Bsize = 0;
    }
    // Fill them
    for (const auto &cp : cell_map)
    {
        const auto &c = cp.second;
        P_MAX = std::max(P_MAX, (int)c.net_list.size());
        for (const auto &netId : c.net_list)
        {
            if (c.partition == 0)
            {
                net_map[netId].Asize++;
            }
            else
            {
                net_map[netId].Bsize++;
            }
        }
    }
    gain_bucket.resize(2 * P_MAX + 1);

    for (auto &cell_pair : cell_map)
    {
        Cell &cell = cell_pair.second;
        if (cell.lock_status == 0)
        { // compute for unlocked cells
            int gain = 0;
            for (const auto &net_id : cell.net_list)
            {
                Net &net = net_map[net_id];
                if (cell.partition == 0)
                { // Cell is in partition A
                    if (net.Asize == 1)
                        gain++;
                    if (net.Bsize == 0)
                        gain--;
                }
                else
                { // Cell is in partition B
                    if (net.Bsize == 1)
                        gain++;
                    if (net.Asize == 0)
                        gain--;
                }
            }
            cell.cell_gain = gain;
            max_gain_index = std::max(max_gain_index, gain + P_MAX);
            min_gain_index = std::min(min_gain_index, gain + P_MAX);
            gain_bucket[gain + P_MAX].push_back(cell.id);
            cell.bucket_iterator = --gain_bucket[gain + P_MAX].end(); // Update bucket iterator
        }
    }
}

void printGainBucketMaxMin()
{
    std::cout << "[INFO] GAIN BUCKET(size=" << gain_bucket.size() << "): ";
    std::cout << "Max Gain Index: " << max_gain_index << ", ";
    std::cout << "Min Gain Index: " << min_gain_index << std::endl;
}

void printGainBucket()
{
    std::cout << "[INFO] Gain Buckets: Under construction " << std::endl;
    std::cout << "[INFO] P_MAX: " << P_MAX << std::endl;
    for (int i = 0; i < gain_bucket.size(); i++)
    {
        std::cout << "Gain Bucket[" << i << "](" << i - P_MAX << "): ";
        for (auto &cellId : gain_bucket[i])
        {
            std::cout << cellId << " ";
        }
        std::cout << std::endl;
    }
}

// ----------------------------------------------------------
// 7) Single FM Pass
// ----------------------------------------------------------
int runOnePass(int mlMode)
{
    int cutBefore = computeCutSize();
    int bestCut = cutBefore;
    int iter = 0;

    // Build the gain buckets
    initializeGainBuckets();
    // printGainBucket();

    // Takeing cell from max to min
    while (max_gain_index >= min_gain_index)
    {
        iter++;
        // printGainBucket();

        int currMaxGain = max_gain_index - P_MAX;
        auto &bucketList = gain_bucket[max_gain_index];

        if (bucketList.empty())
        {
            max_gain_index--;
            continue;
        }

        std::string targetCellId = bucketList.back();
        bucketList.pop_back();
        if (bucketList.empty())
        {
            max_gain_index--;
        }

        auto &targetCell = cell_map[targetCellId];
        if (targetCell.lock_status == 1)
        {
            // skip locked ones
            continue;
        }

        if (mlMode == 1)
        {
            double prob = mlScoreCell(targetCell);
            if (prob < 0.5)
            {
                // lock and skip
                targetCell.lock_status = 1;
                continue;
            }
        }

        // Try moving if it keeps balance
        long long oldA = areaA;
        long long oldB = areaB;

        if (!canMoveCellBalanced(targetCellId))
        {
            // not feasible, lock & revert area
            std::cout << "Caused Imbalance Node:" << targetCellId << std::endl;
            targetCell.lock_status = 1;
            areaA = oldA;
            areaB = oldB;
            continue;
        }

        // Update and move
        updateGainsAfterMove(targetCellId);

        // Move the cell
        targetCell.lock_status = 1;
        int oldPartition = targetCell.partition;
        targetCell.partition = (oldPartition == 0) ? -1 : 0;

        // Update the cut
        int oldCut = cutBefore;
        cutBefore = cutBefore - targetCell.cell_gain;
        std::cout << "NewCut: " << cutBefore << std::endl;
        if (cutBefore < bestCut)
        {
            bestCut = cutBefore;
        }
    }
    return bestCut;
}

// ----------------------------------------------------------
// 8) ML Scoring: A logistic gating that decides if a cell is
//    "likely beneficial" to move
// ----------------------------------------------------------
double mlScoreCell(const Cell &c)
{
    // Features
    double f_area = std::log(1.0 + (double)c.area);
    double f_deg = (double)c.net_list.size();
    double f_gain = (double)c.cell_gain;

    // Hard-coded weights
    double w_area = 0.03;
    double w_deg = 0.05;
    double w_gain = 0.15;
    double bias = -2.0;

    double z = (w_area * f_area + w_deg * f_deg + w_gain * f_gain + bias);
    double prob = 1.0 / (1.0 + std::exp(-z));
    return prob; // 0..1
}

// Helper function to update a cell's gain and move it in the gain bucket
void updateCellGain(Cell &cell, int delta)
{
    int old_gain = cell.cell_gain;
    int new_gain = old_gain + delta;

    // Compute bucket indices
    int old_index = old_gain + P_MAX;
    int new_index = new_gain + P_MAX;

    // Remove from old bucket
    auto &old_bucket = gain_bucket[old_index];
    // printGainBucket();
    old_bucket.erase(cell.bucket_iterator);
    cell.bucket_iterator = gain_bucket[0].end(); // Update bucket iterator

    // Update old bucket becomes empty
    if (old_bucket.empty())
    {
        if (old_index == max_gain_index)
        {
            while (max_gain_index > 0 && gain_bucket[max_gain_index].empty())
            {
                max_gain_index--;
            }
        }
        if (old_index == min_gain_index)
        {
            while (min_gain_index < gain_bucket.size() - 1 && gain_bucket[min_gain_index].empty())
            {
                min_gain_index++;
            }
        }
    }

    // Ensure max_gain_index and min_gain_index are within valid range
    max_gain_index = std::max(0, max_gain_index);
    min_gain_index = std::min((int)gain_bucket.size() - 1, min_gain_index);

    // Add to new bucket
    cell.cell_gain = new_gain;
    auto &new_bucket = gain_bucket[new_index];
    new_bucket.push_back(cell.id); // Update bucket iterator
    cell.bucket_iterator = --new_bucket.end();

    // Update max_gain_index and min_gain_index
    max_gain_index = std::max(max_gain_index, new_index);
    min_gain_index = std::min(min_gain_index, new_index);
}

void updateGainsAfterMove(const std::string &base_cell_id)
{
    Cell &base_cell = cell_map[base_cell_id];
    int from_partition = base_cell.partition;
    int to_partition = from_partition == 0 ? -1 : 0;

    for (const auto &net_id : base_cell.net_list)
    {
        Net &net = net_map[net_id];

        // Check critical nets before the move
        if (from_partition == 0)
        { // Moving from A to B
            if (net.Bsize == 0)
            {
                for (const auto &cell_id : net.cell_list)
                {
                    Cell &cell = cell_map[cell_id];
                    if (cell.lock_status == 0 && cell.id != base_cell_id)
                    {
                        updateCellGain(cell, 1);
                    }
                }
            }
            else if (net.Bsize == 1)
            {
                for (const auto &cell_id : net.cell_list)
                {
                    Cell &cell = cell_map[cell_id];
                    if (cell.lock_status == 0 && cell.partition == 1 && cell.id != base_cell_id)
                    {
                        updateCellGain(cell, -1);
                        break;
                    }
                }
            }
        }
        else
        { // Moving from B to A
            if (net.Asize == 0)
            {
                for (const auto &cell_id : net.cell_list)
                {
                    Cell &cell = cell_map[cell_id];
                    if (cell.lock_status == 0 && cell.id != base_cell_id)
                    {
                        updateCellGain(cell, 1);
                    }
                }
            }
            else if (net.Asize == 1)
            {
                for (const auto &cell_id : net.cell_list)
                {
                    Cell &cell = cell_map[cell_id];
                    if (cell.lock_status == 0 && cell.partition == 0 && cell.id != base_cell_id)
                    {
                        updateCellGain(cell, -1);
                        break;
                    }
                }
            }
        }

        // Update net distribution
        if (from_partition == 0)
        {
            net.Asize--;
            net.Bsize++;
        }
        else
        {
            net.Bsize--;
            net.Asize++;
        }

        // Check critical nets after the move
        if (to_partition == 0)
        { // Moved to A
            if (net.Bsize == 0)
            {
                for (const auto &cell_id : net.cell_list)
                {
                    Cell &cell = cell_map[cell_id];
                    if (cell.lock_status == 0 && cell.id != base_cell_id)
                    {
                        updateCellGain(cell, -1);
                    }
                }
            }
            else if (net.Bsize == 1)
            {
                for (const auto &cell_id : net.cell_list)
                {
                    Cell &cell = cell_map[cell_id];
                    if (cell.lock_status == 0 && cell.partition == 1 && cell.id != base_cell_id)
                    {
                        updateCellGain(cell, 1);
                        break;
                    }
                }
            }
        }
        else
        { // Moved to B
            if (net.Asize == 0)
            {
                for (const auto &cell_id : net.cell_list)
                {
                    Cell &cell = cell_map[cell_id];
                    if (cell.lock_status == 0 && cell.id != base_cell_id)
                    {
                        updateCellGain(cell, -1);
                    }
                }
            }
            else if (net.Asize == 1)
            {
                for (const auto &cell_id : net.cell_list)
                {
                    Cell &cell = cell_map[cell_id];
                    if (cell.lock_status == 0 && cell.partition == 0 && cell.id != base_cell_id)
                    {
                        updateCellGain(cell, 1);
                        break;
                    }
                }
            }
        }
    }
}

void printPartition(const std::string &partitionName, bool showLockStatus = false)
{
    std::cout << "Partition " << partitionName << ": ";
    for (auto &cp : cell_map)
    {
        if (cp.second.partition == 0 && partitionName == "A")
        {
            if (showLockStatus)
                std::cout << "{ " << cp.first << "," << cp.second.lock_status << "} ";
            else
                std::cout << cp.first << " ";
        }

        if (cp.second.partition == -1 && partitionName == "B")
        {
            if (showLockStatus)
                std::cout << "{ " << cp.first << "," << cp.second.lock_status << "} ";
            else
                std::cout << cp.first << " ";
        }
    }
    std::cout << std::endl;
}

// main()
int main(int argc, char **argv)
{
    if (argc < 5)
    {
        std::cerr << "Usage: " << argv[0]
                  << " <benchmarkName> <areaMode=1|2> <mlMode=0|1> <maxPasses>\n";
        std::cerr << "Example: " << argv[0] << " superblue18 2 1 5\n";
        return 1;
    }

    std::string benchmarkName = argv[1];
    int areaMode = std::stoi(argv[2]);  // 1 => all gates area=1, 2 => width*height
    int mlMode = std::stoi(argv[3]);    // 0 => ML off, 1 => ML on
    int maxPasses = std::stoi(argv[4]); // number of passes

    std::string nodeFile = "data/" + benchmarkName + "/" + benchmarkName + ".nodes";
    std::string netFile = "data/" + benchmarkName + "/" + benchmarkName + ".nets";

    auto startTP = std::chrono::system_clock::now();
    std::time_t startT = std::chrono::system_clock::to_time_t(startTP);
    std::cout << "[INFO] Starting Partition at " << std::ctime(&startT) << std::endl;

    // 1) Read Cells
    readCells(nodeFile, areaMode);

    // 2) Random Initial Partition
    initPartition();

    // 3) Read Nets
    readNets(netFile);

    // 4) Compute initial cut
    int currentCut = computeCutSize();
    bestCutGlobal = currentCut;
    bestAreaA = areaA;
    bestAreaB = areaB;

    // Save initial partition as best
    bestPartitionRecord.clear();
    for (const auto &cp : cell_map)
    {
        bestPartitionRecord[cp.first] = cp.second.partition;
    }

    std::cout << "[INFO] Initial Cut Size = " << currentCut << std::endl;

    // 5) Multi-pass FM
    for (int pass = 1; pass <= maxPasses; pass++)
    {
        // Unlock all cells
        for (auto &cp : cell_map)
        {
            cp.second.lock_status = 0;
        }
        // std::cout << " [PASS " << pass << "] Starting Pass..." << std::endl;
        int passCut = runOnePass(mlMode);

        if (passCut < bestCutGlobal)
        {
            bestCutGlobal = passCut;
            // record new best partition
            bestAreaA = areaA;
            bestAreaB = areaB;
            for (const auto &cp : cell_map)
            {
                bestPartitionRecord[cp.first] = cp.second.partition;
            }
            std::cout << " [PASS " << pass << "] Improved Cut => " << passCut << std::endl;
        }
        else
        {
            // no improvement => revert to best partition and break
            std::cout << " [PASS " << pass << "] No improvement (cut=" << passCut
                      << "), revert & stop.\n";
            // revert
            areaA = bestAreaA;
            areaB = bestAreaB;
            for (auto &cp : cell_map)
            {
                cp.second.partition = bestPartitionRecord[cp.first];
                cp.second.lock_status = 0;
            }
            break;
        }
    }

    // Print Partition A:
    // printPartition("A");

    // Print Partition B:
    // printPartition("B");

    // 6) Done. bestCutGlobal is the final cut.
    std::cout << "[INFO] Final Cut Size   = " << bestCutGlobal << std::endl;
    std::cout << "[INFO] Final partition areas: A=" << bestAreaA
              << ", B=" << bestAreaB << std::endl;

    auto endTP = std::chrono::system_clock::now();
    std::time_t endT = std::chrono::system_clock::to_time_t(endTP);
    double elapsedSec = std::chrono::duration<double>(endTP - startTP).count();

    std::cout << "[INFO] Ending Partition at " << std::ctime(&endT) << std::endl;
    std::cout << "[INFO] Elapsed time: " << elapsedSec << " seconds" << std::endl;

    // Counting cells in A/B
    long countA = 0, countB = 0;
    for (const auto &cp : bestPartitionRecord)
    {
        if (cp.second == 0)
            countA++;
        else
            countB++;
    }
    std::cout << "[INFO] #Cells in A=" << countA
              << "  #Cells in B=" << countB << std::endl;

    return 0;
}