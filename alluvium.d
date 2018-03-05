// **********************************************************************************************
//
// alluvium.d - a sediment accumulation model for an alluvial fan
//
// by Walt McNab
//
// Key conceptual components:
// (1) Rivulet generation scheme in response to land surface gradient
// (2) Representation of physical-empirical sediment deposition 
//
// **********************************************************************************************

// libraries
import std.stdio;           // I/O and file system
import std.math;            // math function library
import std.algorithm;       // algorithm for working with ranges, sequences, etc.
import std.array;           // array operations
import std.string;          // string function support
import std.conv;            // automatic conversions between types
import std.typecons;        // misc. meta-programming support, including type constructors, tuples, etc.
import std.range;           // enables range functions
import std.mathspecial;     // for inverse normal distribution function
import std.random;          // for random number generation
import utility;             // special defined functions, co-located with alluvium.d
    
// physical constants
const(double) g = 9.807;                    // gravitational acceleration (m/sec^2)
const(double) rhoW = 1000;                  // density of water (kg/m^3)
const(double) rhoS = 2650;                  // density of sedimentary particles (kg/m^3)
const(double) rhob = 1600;                  // sediment bulk density (kg/m^3)
const(double) s = 2.65;                     // sediment specific gravity
const(double) uW = 8.9e-4;                  // dynamic viscosity of water (Pa*sec)
const(double) nu = 8.9e-7;                  // kinematic viscosity (m^2/sec)    
    
    
// **********************************************************************************************
//
// sedimentology classes & supporting functions
//
// Key assumptions:
// (1) Water column in each cell in stream is perfectly stirred by turbulent mixing
// (2) Flow, and hence sediment deposition rates, jumps from one steady-state condition to next
// (3) Water is conserved in stream (no percolation loss)
//
// **********************************************************************************************   


class Drainage{

    // characteristics of a drainage contact area (location, influx, sediment composition and variability)
    int numSeds, disCellIndex, edgeSrcIndex;
    double fracQ, fracProb, fracMult;
    double[] fracMin;
    double[] fracMax;
    double[] bandwidth;
    
    // constructor
    this(double xs, double ys, double fracQ, double fracProb, double fracMult, double[] fracMin, double[] fracMax, Grid grid){
        this.numSeds = fracMin.length;                  // infer number of sediment classes
        auto src = grid.DischargeCell(xs, ys);
        this.disCellIndex = src[0];                     // grid index of discharge source
        this.edgeSrcIndex = src[1];                     // edge index of discharge source   
        this.fracQ = fracQ;                             // fraction of reference discharge passing through this drainage
        this.fracProb = fracProb;                       // probability parameter for sediment composition random walk   
        this.fracMult = fracMult;                       // multiplier parameter for sediment composition random walk        
        this.fracMin = fracMin;                         // constraints on mineralogical composition
        this.fracMax = fracMax;
        bandwidth.length = numSeds;
        bandwidth[] = fracMax[] - fracMin[];
        this.bandwidth = bandwidth;
        }

    bool CheckValid(double[] propSet){
        // check that a proposed composition is within proportion bounds
        bool valid = true;
        for (int i = 0; i < numSeds; ++i){
            if ((propSet[i] < fracMin[i]) || (propSet[i] > fracMax[i])){
                valid = false;
                break;
                }
            }
        return valid;
        }

    Tuple!(double[], int) Seed(){
        // random sediment proportion makeup for starting point
        bool valid = false;
        int perturbSelect;
        double[] seed;
        seed.length = numSeds;
        while (valid != true){
            double[] raw;
            for (int i = 0; i < numSeds; ++i){raw ~= uniform01() * bandwidth[i];}
            seed[] = Normalize(raw);
            valid = CheckValid(seed);
            }
        perturbSelect = to!int(floor(uniform01()*numSeds));     // sediment class selection seed
        return tuple(seed, perturbSelect);
        }

    Tuple!(double[], int) RandomWalk(int perturbSelect, double[] propSet){
        // incremental random walk in proportion space
        double rSelect;
        double[] adjSet;
        double[] adjSetNorm;
        double[] perturb;
        bool valid = false;
        adjSet.length = numSeds;
        perturb.length = numSeds;
        while (valid != true){
            rSelect = uniform01();
            if (rSelect >= fracProb){perturbSelect = to!int(floor(uniform01()*numSeds));}
            perturb[] = 0.;
            double[] r;
            for (int i = 0; i < numSeds; ++i){r ~= uniform01() * bandwidth[i] * fracMult;}
            perturb[perturbSelect] += r[perturbSelect]; 
            adjSet[] = propSet[] + perturb[];
            adjSetNorm = Normalize(adjSet);
            valid = CheckValid(adjSetNorm);
            }
        return tuple(adjSetNorm, perturbSelect);       
        }
        
    }

class Discharge{

    // discharge (stream flow) intensity and duration

    double logMeanQ, logStdevQ, debrisFlowQ, dryCutoffQ, a, b, logtRes, endTime;

    // constructor
    this(){
        string lineInput;
        string[][] parsedInput; 
        auto inputFile = File("discharge.txt","r");             // read in contents of input file
        while (!inputFile.eof()) {
            lineInput = inputFile.readln();
            parsedInput ~= split(lineInput);
            }
        inputFile.close();          
        this.logMeanQ = to!double(parsedInput[0][1]);           // log mean discharge (log units)
        this.logStdevQ = to!double(parsedInput[1][1]);          // log standard deviation discharge (log units)
        this.debrisFlowQ = to!double(parsedInput[2][1]);        // cut-off maximum discharge, above which debris flows occur
        this.dryCutoffQ = to!double(parsedInput[3][1]);         // cut-off minimum discharge, to prevent prolonged low-flow
        this.a = to!double(parsedInput[4][1]);                  // log dt axis intercept
        this.b = to!double(parsedInput[5][1]);                  // slope of the log dt vs log Q relationship
        this.logtRes = to!double(parsedInput[6][1]);            // log standard deviation of event duration, to provide "wiggle"
        this.endTime = to!double(parsedInput[7][1]);            // end-of-simulation time
        writeln("Read hydrologic generator parameters.");       
        }

    Tuple!(double, double) Event(){
        // generate a random discharge event (magnitude and duration)
        double logQ, dt, logdt, p1, p2;
        double Q = 1e-10;
        while(Q < dryCutoffQ){
            p1 = uniform01();
            p2 = uniform01();       
            logQ = logMeanQ + logStdevQ*normalDistributionInverse(p1);
            logdt = a + b*logQ + logtRes*normalDistributionInverse(p2);
            Q = 10.^^logQ;
            dt = 10.^^logdt;
            }
        return tuple(Q, dt);
        }
    
    }
    
    
class Sediment {

    // properties and methods (basic physics) associated with different sediment sizes

    string name;
    double d50, d90, D, vs;

    // constructor
    this(string name, double d50, double d90){
        this.name = name;                   // name of sediment classes
        this.d50 = d50;                     // diameter of 50th percentile
        this.d90 = d90;                     // diameter of 90th percentile      
        this.D = dnd();                     // dimensionless particle diameter
        this.vs = SettlingVel();            // settling velocity        
        }
        
    double dnd(){
        // non-dimensional particle diameter
        return ((s-1.)*g/nu^^2)^^(1./3.) * d50;
        }

    double SettlingVel(){
        // particle settling velocity (Cheng, 1997 approximation)
        return nu/d50 * ((25. + 1.2*D^^2)^^0.5 - 5.)^^1.5;
        }       
        
    }


class Stream {

    // properties and methods (e.g., flow hydraulics, sediment load, deposition rates) for a stream as a whole

    double w, dx, A, n, xBL, dsMax, S0;
    
    // constructor
    this(double w){
        string lineInput;
        string[][] parsedInput; 
        auto inputFile = File("stream.txt","r");                // read in contents of input file
        while (!inputFile.eof()) {
            lineInput = inputFile.readln();
            parsedInput ~= split(lineInput);
            }
        inputFile.close();
        this.w = w;                                         // stream width
        this.dx = w;                                        // cell (a.k.a.reach); equal to "w" for 2-D grid
        this.A = w*dx;                                      // areal footprint of a cell
        this.n = to!double(parsedInput[0][1]);              // Manning roughness coefficient (global)
        this.xBL = to!double(parsedInput[1][1]);            // bed-load layer thickness maximum
        this.dsMax = to!double(parsedInput[2][1]);          // maximum implied change in dSed/h before adjusting flow rate/new time step
        this.S0 = to!double(parsedInput[3][1]);             // slope of stream as it exits upland drainage area and enters alluvium
        writeln("Read stream parameters.");     
        }
    
    double[][] DebrisFlow(Cell[] reach, double Q, double[] drainageFrac){
        // calculate debris flow thickness (~ flow thickness --> ad hoc assumption)
        static int iterMax = 3;                     // fixed number of iterations used to approximate Manning equation
        double[][] J;
        double[] S;
        double[] h;
        int numCells = reach.length;
        h.length = numCells;
        J.length = numCells;
        h[] = 0.;       
        Tuple!(double, double, double) flow;
        for (int iter = 0; iter < iterMax; ++iter){         // iteratively approximate Manning equation solution
            S = Slope(numCells, h, reach);                  // slope = slope of water surface
            for (int i = 0; i < numCells; ++i){
                flow = Manning(Q, max(S[i], S0));
                h[i] = flow[0];
                }
            }   
        for (int i = 0; i < numCells; ++i){
            for (int j = 0; j < drainageFrac.length; ++j){
                J[i] ~= drainageFrac[j] * h[i];
                }
            }
        return J;
        }
        
    Tuple!(double[][], double) SettlingEvent(Cell[] reach, double Q, Sediment[] seds, double[] drainageFrac){

        // manager function for stream deposition event; returns thickness of each sediment size class added to each cell/reach
        double h0, hBL0r, q0, C0SL, C0BL, CcSL, CcBL, jSL, jBL;
        double[] S;
        double[] h;
        double[] hBLr;
        double[] q;
        double[][] J;
        double dtImply = 1e+30;                     // time step size required to scale maximum deposition rate to constraint
        static double epsilon = 1e-30;              // factor used to prevent divide-by-zero when computing time step size
        static int iterMax = 3;                     // fixed number of iterations used to approximate Manning equation
        
        int numCells = reach.length;                // initialize variables used with Manning flow equation
        h.length = numCells;                    
        hBLr.length = numCells;
        q.length = numCells;
        h[] = 0.;
        hBLr[] = 0.;
        q[] = 0.;
        
        Tuple!(double, double) C;                   // container variables to handle various function outputs
        Tuple!(double, double) FS;      
        Tuple!(double, double) FB;      
        Tuple!(double, double, double) flow;
        J.length = numCells;
        
        // *** step #1: calculate flow velocities, based on conditions at start of time step ***
        auto flow0 = Manning(Q, S0);                // incoming flow from uplands
        h0 = flow0[0];                              // water depth
        hBL0r = flow0[1]/h0;                        // bed load boundary depth/water depth ratio
        q0 = flow0[2];                              // stream flow velocity

        // iteratively approximate Manning equation solution
        for (int iter = 0; iter < iterMax; ++iter){
            S = Slope(numCells, h, reach);              // slope = slope of water surface
            for (int i = 0; i < numCells; ++i){
                flow = Manning(Q, S[i]);
                h[i] = flow[0];
                hBLr[i] = flow[1]/flow[0];
                q[i] = flow[2];
                }
            }
            
        for(int j = 0; j < seds.length; ++j){       // for each size class ...
        
            // step #2a: calculate the equilibrium sediment concentrations (incoming)
            C = Capacity(Q, drainageFrac[j], h0, hBL0r, q0, seds[j]);
            C0SL = C[0];
            C0BL = C[1];
            
            for (int i = 0; i < numCells; ++i){
            
                // step #2b: calculate the equilibrium sediment concentrations (top layer in each cell/reach)           
                C = Capacity(Q, reach[i].lithFrac[0][j], h[i], hBLr[i], q[i], seds[j]);
                CcSL = C[0];
                CcBL = C[1];
                
                // step #3: calculate the total deposition rate (per reach, per sediment size class)
                FS = DepRate(C0SL, Q, h[i], q[i], CcSL, seds[j]);                   // suspended load component
                jSL = FS[0];
                C0SL = FS[1];
                
                FB = DepRate(C0BL, Q*hBLr[i], h[i], q[i], CcBL, seds[j]);   // bed load component
                jBL = FB[0];
                C0BL = FB[1];
                J[i] ~= (jSL + jBL) / rhob;     // note that units are length/time
                
                }
            }
            
        // implied time step to reach specified thickness
        for (int i = 0; i < numCells; ++i){dtImply = min(dsMax*h[i]/(sum(J[i])+epsilon), dtImply);}
        return tuple(J, dtImply);  
        }   
        
    Tuple!(double, double) Capacity(double Q, double f, double h, double hBRr, double q, Sediment sed){
        // calculate the "equilibrium" sediment concentration(s) for the given particle size for given cell/reach
        double cSL = f * SuspendedLoad(q, h, sed) * w / Q;
        double cBL = f * BedLoad(q, h, sed) * w / (Q*hBRr);     
        return tuple(cSL, cBL);
        }   
        
    Tuple!(double, double) DepRate(double C0, double Qf, double h, double q, double Cc, Sediment sed){
        // calculate deposition rate across cells in stream, for suspended OR bed load, for a particular size class
        // note that Qf is the effective flow rate for suspended load (fill stream depth) or bed load (partial layer)
        double Cdep, Ceff;                      // sediment mass concentrations: available for deposition, and effective
        double J;                               // deposition rate (mass/time)
        Cdep = max(C0-Cc, 0.);                  // calculate difference in carrying capacity; a drop means deposition (erosion not modeled)
        Ceff = Cdep * Qf/(sed.vs*A + Qf);       // deposition rate reflects "surplus" sediment concentration, settling velocity, areal footprint, and discharge
        J = sed.vs*Ceff;                        // deposition rate expressed as mass area^-1 time^-1
        C0 = min(Cc+Ceff, C0);                  // outflow conc is equal to capacity conc plus effective conc (deposition), or inflow conc (no deposition)
        return tuple(J, C0);
        }
        
    double Me(double q, double ucr, Sediment sed){
        // mobility parameter for simplified sediment load models (for all class sizes)
        double me;
        if (q > ucr){me = (q-ucr)/((s-1)*g*sed.d50)^^0.5;}
        else {me = 0;}
        return me;
        }       

    double Ucr(double h, Sediment sed){
        // critical velocity for currents (for all class sizes)
        double ucr;
        if (sed.d50 < 0.0005){ucr = 0.19*sed.d50^^0.1 * log10(12*h/(3.*sed.d90));}
        else {ucr = 8.5*sed.d50^^0.6 * log10(12*h/(3.*sed.d90));}
        return ucr;
        }
            
    double BedLoad(double q, double h, Sediment sed){
        // simplified bed load transport formula for steady-state current flow (for all class sizes)
        double ucr = Ucr(h, sed);
        double me = Me(q, ucr, sed);
        return 0.015 * rhoS * q * h * (sed.d50/h)^^1.2 * me^^1.5;   // dimensions = mass length^-1 time^-1
        }

    double SuspendedLoad(double q, double h, Sediment sed){
        // simplified suspended load transport formula for steady-state current flow
        double ucr = Ucr(h, sed);
        double me = Me(q, ucr, sed);                            
        return 0.012 * rhoS * q * sed.d50 * me^^2.4 * sed.D^^(-0.6);    // dimensions = mass length^-1 time^-1
        }       
        
    double[] Slope(int numCells, double[] h, Cell[] reach){
        // calculate slopes along the stream (approximated by streambed elevation slope, as opposed to water surface)
        double[] S;
        for (int i = 0; i < numCells-1; ++i){
            S ~= (reach[i].zTop + h[i] - reach[i+1].zTop - h[i+1])/dx;
            }
        S ~= (reach[numCells-1].zTop + h[numCells-1]) / dx;     // last cell
        return S;
        }
        
    Tuple!(double, double, double) Manning(double Q, double S){
        // water depth and flow velocity at a given location via the Manning equation
        double h, hBL, q;
        h = 10.^^(3./5. * log10(Q/w * n * 1./sqrt(S)));             // stream depth
        hBL = min(h, xBL);                                          // bed-load layer thickness
        q = Q/(w*h);                                                // flow velocity
        return tuple(h, hBL, q);
        }   
    
    }


// *** input functions ***


Drainage[] ReadDrainage(Grid grid){
    // read in drainage sediment makeup
    int lineNum = 0;
    double xs, ys, fracQ, fracProb, fracMult;
    double[] fracMin;
    double[] fracMax;
    string lineInput;
    string[][] parsedInput; 
    Drainage[] drainageObj;
    auto inputFile = File("drainage.txt","r");
    while (!inputFile.eof()) {
        lineInput = inputFile.readln();
        parsedInput ~= split(lineInput);
        }
    inputFile.close();      
    fracMin.length = parsedInput[1].length - 7;             // number of sediment size classes
    fracMax.length = parsedInput[1].length - 7; 
    for (int i = 1; i < parsedInput.length; ++i){
        if (parsedInput[i][6]=="min"){
            xs = to!double(parsedInput[i][1]);
            ys = to!double(parsedInput[i][2]);
            fracQ = to!double(parsedInput[i][3]);
            fracProb = to!double(parsedInput[i][4]);
            fracMult = to!double(parsedInput[i][5]);            
            for (int j = 7; j < parsedInput[i].length; ++j){
                fracMin[j-7] = to!double(parsedInput[i][j]);
                }
            }
        else{
            for (int j = 7; j < parsedInput[i].length; ++j){
                fracMax[j-7] = to!double(parsedInput[i][j]);
                }   
                drainageObj ~= new Drainage(xs, ys, fracQ, fracProb, fracMult, fracMin, fracMax, grid);
            }
        }
    writeln("Read drainage physical characteristics and sediment composition constraints.");    
    return drainageObj;
    }

Sediment[] ReadSeds(){
    // read in sediment class characteristics
    int lineNum = 0;
    double d50, d90;
    string lineInput, name;
    string[] parsedInput;
    Sediment[] seds;
    auto inputFile = File("sediments.txt","r");
    while (!inputFile.eof()) {
        lineInput = inputFile.readln();
        if (lineNum > 0) {                                  // ignore header line in input file 
            parsedInput = split(lineInput);
            name = parsedInput[0];
            d50 = to!double(parsedInput[1]);            
            d90 = to!double(parsedInput[2]);            
            Sediment seds_member = new Sediment(name, d50, d90);
            seds ~= seds_member;
            }
        ++lineNum;
        }
    inputFile.close();  
    writeln("Read sediment properties.");
    return seds;
    }           
    
    
// **********************************
//
// grid-oriented class definitions 
//
// **********************************   
    

class Grid{
    
    int nx, ny;
    double d, f, dz, dzTop;
    int[4] dirOpp;
    Cell[] cell;
    
    this(int[] col, int[] row, double[] z, int numSeds){
        // read grid property values and cell elevations
        int icol, jrow;
        double x0, y0, x, y;
        double[] lithFrac0;             
        string lineInput;
        string[][] parsedInput; 
        auto inputFile = File("grid_params.txt","r");           // read in contents of input file
        while (!inputFile.eof()) {
            lineInput = inputFile.readln();
            parsedInput ~= split(lineInput);
            }
        inputFile.close();          
        x0 = to!double(parsedInput[0][1]);              // grid origin
        y0 = to!double(parsedInput[1][1]);      
        this.d = to!double(parsedInput[2][1]);          // grid spacing (must be the same for both axes)
        this.f = to!double(parsedInput[3][1]);          // incoming flow bias multiplier for choosing gradient  
        dz = to!double(parsedInput[4][1]);
        dzTop = to!double(parsedInput[5][1]);
        for (int i = 1; i <= numSeds; ++i){             // default base sediment composition
            lithFrac0 ~= to!double(parsedInput[6][i]);
            }       
        this.nx = MaxIntArray(col) + 1;                 // number of cells along each axis
        this.ny = MaxIntArray(row) + 1;        
        this.dirOpp = [2, 3, 0, 1];                     // opposing directions (for incoming flow bias)
        for (int i = 0; i < col.length; ++i){           // * create individual cell objects within grid *
            icol = col[i];
            jrow = row[i];
            x = x0 + (icol+0.5)*d;
            y = y0 + (jrow+0.5)*d;
            this.cell ~= new Cell(x, y, z[i], z[i]+dzTop, icol, jrow, nx, ny, dz, dzTop, 1, lithFrac0);
            }
        writeln("Populated grid object.");
        }

    Tuple!(int, int) DischargeCell(double xs, double ys){
        // find cell index number corresponding to discharge point (xs, ys); used by other classes
        int jrow, icol, disCellIndex, edgeSrcIndex;
        bool[4] edgeContainer;
        jrow = to!int(floor(ys/d));
        icol = to!int(floor(xs/d));
        disCellIndex = jrow * nx + icol;
        edgeContainer = [jrow==ny-1, icol==nx-1, jrow==0, icol==0];
        edgeSrcIndex = ArrayIndex(edgeContainer, true);      // edge containing discharge point
        return tuple(disCellIndex, edgeSrcIndex);
        }

    double[] GetGrad(int cellIndex){
        // gradient is positive for flow directed out of cell
        double dzdx;
        double[] grad;
        foreach (cn; cell[cellIndex].connects[1]){
            dzdx = (cell[cellIndex].zTop - cell[cn].zTop) / d;
            grad ~= dzdx;
            }
        return grad;
        }
        
    bool CheckEdge(int cellIndex, int edgeSrcIndex){
        // check to see if cell has reached a non-discharging model edge
        bool edgeState = false;       // default assumption
        if (cell[cellIndex].edgeIndex >= 0){         // it's along an edge
            if (cell[cellIndex].edgeIndex != edgeSrcIndex){ 
                // it's not the same edge as the discharge source
                edgeState = true;
                }
            }
        return edgeState;
        }       
        
    void FillHole(int cellIndex){
        // fill in any depressions in ground surface found by March method with finest-grained material
        int adjIndex;
        double[] adjElev;
        double[] dSed;
        dSed.length = cell[cellIndex].lithFrac[0].length;       // number of sediment size classes (since not passed to method)
        dSed[] = 0.;
        for (int i = 0; i < cell[cellIndex].connects[1].length; ++i){
            adjIndex = cell[cellIndex].connects[1][i];
            adjElev ~= cell[adjIndex].zTop;
            }
        double newElev = sum(adjElev)/adjElev.length;           // corrected elevation = average of those of surrounding cells
        dSed[0] = newElev - cell[cellIndex].zTop;
        cell[cellIndex].Deposit(dSed);                          // add/mix material to cell
        writeln("\tFilled in low point in grid.");
        }   
        
    Tuple!(int, int, bool) March(int cellIndex, int p, double dirF){
        // pick a new cell in path, based on gradient and momentum considerations
        int flowThruIndex, flowSelect, flowDir, newCell;
        bool valid = true;                              // assume addition pathway is valid (no holes) until proven otherwise
        double[] grad;
        double[] outflows;
        grad = GetGrad(cellIndex);
        for (int i = 0; i < grad.length; ++i){          // find potential positive fluxes out of cell
            if (grad[i]>=0.){outflows ~= grad[i];}
            else {outflows ~= 0.;}
            }
        // check if path gets stuck in a local minimum  
        if (sum(outflows)==0.){
            FillHole(cellIndex);        // fill in depression
            valid = false;
            }   
        flowThruIndex = ArrayIndex(cell[cellIndex].connects[0], dirOpp[p]);     // determine if flow-through possible
        if (flowThruIndex > -1){outflows[flowThruIndex] *= dirF;}               // add optional weight to flow through direction
        flowSelect = SelectRand(outflows);                                      // pick an index at random for new flow direction
        flowDir = cell[cellIndex].connects[0][flowSelect];                      // outgoing flow direction
        newCell = cell[cellIndex].connects[1][flowSelect];                      // outgoing flow cell
        return tuple(flowDir, newCell, valid);      
        }

    int[] Path(int disCellIndex, int edgeSrcIndex){
        // compute random pathway from discharge point to edge
        int flowDir, newCellIndex;
        int[] path;
        Tuple!(int, int, bool) flow;
        path ~= disCellIndex;                                       // every path starts at the discharge location cell
        flow = March(disCellIndex, 0, 1);                           // initiate cell-by-cell migration
        flowDir = flow[0];
        newCellIndex = flow[1];
        path ~= newCellIndex;
        while (CheckEdge(newCellIndex, edgeSrcIndex) == false){
            flow = March(newCellIndex, flowDir, f);
            if (flow[2]==true){                                     // only append to path if March reports valid addition
                flowDir = flow[0];
                newCellIndex = flow[1];
                path ~= newCellIndex;
                }
            }
        return path;
        }

    void WriteSedOuput(Sediment[] seds){
        // write sediment stack summary to tabular text output file
        double z;
        string header, outputLine;
        // open file and write header
        auto outputFile = File("sed_distribution.csv", "w");    
        header = "x" ~ "," ~ "y" ~ "," ~ "z";
        foreach(sed; seds){header ~= "," ~ sed.name;}
        outputFile.writeln(header); 
        for (int i = 0; i < cell.length; ++i){
            for (int j = 0; j < cell[i].numLayers-1; ++j){      // for all layers beneath surface layer
                z = cell[i].z0 + (j+0.5)*cell[i].dz;
                outputLine = to!string(cell[i].x) ~ "," ~ to!string(cell[i].y) ~ "," ~ to!string(z);
                for (int k = 0; k < seds.length; ++k){outputLine ~= "," ~ to!string(cell[i].lithFrac[j][k]);}
                outputFile.writeln(outputLine);
                }
            // for surface layer                
            z = cell[i].zTop - 0.5*cell[i].dzTop;   
            outputLine = to!string(cell[i].x) ~ "," ~ to!string(cell[i].y) ~ "," ~ to!string(z);
            for (int k = 0; k < seds.length; ++k){outputLine ~= "," ~ to!string(cell[i].lithFrac[cell[i].numLayers-1][k]);}
            outputFile.writeln(outputLine);
            }
        outputFile.close(); 
        writeln("Wrote sediments output file.");
        }   
    
    void WriteElevOuput(Sediment[] seds){
        // output file for land surface elevation and composition distribution
        string header, outputLine;
        auto outputFile = File("land_surface.csv", "w");    
        header = "x" ~ "," ~ "y" ~ "," ~ "z";
        foreach(sed; seds){header ~= "," ~ sed.name;}
        outputFile.writeln(header); 
        for (int i = 0; i < cell.length; ++i){
            outputLine = to!string(cell[i].x) ~ "," ~ to!string(cell[i].y) ~ "," ~ to!string(cell[i].zTop);
            for (int k = 0; k < seds.length; ++k){outputLine ~= "," ~ to!string(cell[i].lithFrac[cell[i].numLayers-1][k]);}
            outputFile.writeln(outputLine);
            }
        outputFile.close(); 
        writeln("Wrote land surface elevation output file.");       
        }
    
    void WriteSurface(){
        // output file for land surface elevation only (for debugging; normally not called)
        string header, outputLine;
        auto outputFile = File("check_surface.csv", "w");   
        header = "x" ~ "," ~ "y" ~ "," ~ "z";
        outputFile.writeln(header); 
        for (int i = 0; i < cell.length; ++i){
            outputLine = to!string(cell[i].x) ~ "," ~ to!string(cell[i].y) ~ "," ~ to!string(cell[i].zTop);
            outputFile.writeln(outputLine);
            }
        outputFile.close(); 
        writeln("Wrote land surface elevations to check file.");        
        }
    
    }

    
class Cell {

    // location, sediment pile (array), and methods for a given cell

    double x, y, z0, zTop, dz, dzTop;
    int numLayers, edgeIndex;
    double[] lithFracInitial;
    double[][] lithFrac;    
    bool[4] edgeContainer;
    Tuple!(int[], int[]) connects;
    
    // constructor
    this(double x, double y, double z0, double zTop, int icol, int jrow, int numCols, int numRows, 
        double dz, double dzTop, int numLayers, double[] lithFrac0){
        this.x = x;                                     // center-of-cell location
        this.y = y; 
        this.z0 = z0;                                   // base (i.e., bedrock) elevation; fixed quantity
        this.zTop = zTop;                               // elevation of top of sediment pile; dynamic quantity
        this.dz = dz;                                   // full slice thickness
        this.dzTop = dzTop;                             // thickness of top layer, which may be incomplete
        this.numLayers = numLayers;                     // current number of sediment layers (full + incomplete top)
        for (int i = 0; i < lithFrac0.length; ++i){     // array of volume proportions of different sediment classes represented in each slice
            lithFracInitial ~= lithFrac0[i];
            }       
        this.lithFrac ~= lithFracInitial;   
        this.connects = AssignConnects(icol, jrow, numCols, numRows);
        edgeContainer = [jrow==numRows-1, icol==numCols-1, jrow==0, icol==0];           // check if cell is along an edge; mark which one 
        this.edgeIndex = ArrayIndex(edgeContainer, true);
        }

    Tuple!(int[], int[]) AssignConnects(int icol, int jrow, int numCols, int numRows){
        // assign connection information: directions and cell index numbers
        int[] connectDir;
        int[] connectIndex;
        if (jrow != numRows-1){
            connectDir ~= 0;
            connectIndex ~= GridIndex(jrow+1, icol, numCols);
            }
        if (icol != numCols-1){
            connectDir ~= 1;
            connectIndex ~= GridIndex(jrow, icol+1, numCols);       
            }
        if (jrow != 0){
            connectDir ~= 2;
            connectIndex ~= GridIndex(jrow-1, icol, numCols);       
            }
        if (icol != 0){
            connectDir ~= 3;
            connectIndex ~= GridIndex(jrow, icol-1, numCols);       
            }
        return tuple(connectDir, connectIndex);
        }

    int GridIndex(int jrow, int icol, int numCols){
        // fetch cell index number corresponding to jrow and icolumn
        return jrow * numCols + icol;
        }
        
    void Deposit(double[] dSed){
        // local deposition of sediments onto stack; re-discretize as warranted
        double layerAugment, extraAdded, newTopLayer, mixFrac;
        double[] dSedComp;
        int numNewFullLayers;
        double totalAdded = sum(dSed);              // sum of new sediment deposit thickness (across all grain size classes)
        dSedComp = Normalize(dSed);                 // normalized particle size distribution of new deposit
        if (dzTop + totalAdded <= dz){
            // add all of the new accumulation to existing top layer
            mixFrac = dzTop/(totalAdded + dzTop);
            lithFrac[numLayers-1] = Blend(mixFrac, lithFrac[numLayers-1], dSedComp);
            dzTop += totalAdded;
            zTop += totalAdded;
            }
        else{
            mixFrac = dzTop/dz;
            lithFrac[numLayers-1] = Blend(mixFrac, lithFrac[numLayers-1], dSedComp);
            layerAugment = dz - dzTop;                                              // amount consumed to fill out old top layer to max thickness           
            extraAdded = totalAdded - layerAugment;                                 // amount leftover after top layer fill-out
            numNewFullLayers = to!int(floor(extraAdded/dz));                        // number of full-thickness layers, above old top
            dzTop = extraAdded - numNewFullLayers*dz;                               // update incomplete top layer thickness (new)
            numLayers += numNewFullLayers + 1;                                      // update the total number of layers
            zTop = z0 + (numLayers-1)*dz + dzTop;                                   // update top-of-sediment elevation 
            for (int i = 0; i < numNewFullLayers+1; ++i){lithFrac ~= dSedComp;}     // add lithology mix to new layers
            }   
        }
    
    }   

    
// *** input function ***


Tuple!(int[], int[], double[]) ReadCells(){
    // read cell elevations from input file
    int[] col;
    int[] row;
    double[] z;
    int lineNum = 0;
    string lineInput;
    string[] parsedInput;   
    auto inputFile = File("cells.txt","r");
    while (!inputFile.eof()){
        lineInput = inputFile.readln();
        if (lineNum > 0) {                      // header
            parsedInput = split(lineInput);
            col ~= to!int(parsedInput[0]);
            row ~= to!int(parsedInput[1]);            
            z ~= to!double(parsedInput[2]);
            }           
        ++lineNum;
        }
    inputFile.close();
    writeln("Read cell base elevations.");
    return tuple(col, row, z);
    }
    
    
// ********
//
//   main 
//
// ********


double FixTime(double dtProp, double t, double tEnd){
    // correct (variable) time step size to prevent implied time exceedances
    double dt;
    if ((t+dtProp) > tEnd){dt = tEnd-t;}        // dt > remaining time, so let dt = remaining time
    else {dt = dtProp;}                         // still room for dt under remaining time; no change
    return dt;
    }   
    
void main(){

    double tStep, dtStep, dtEvent, Q, Qeff;
    double t = 0.;
    int[] sedSelect;    
    double[] dSed;
    double[][] J;
    double[][] drainageFrac;
    string dischargeLine;
    Tuple!(double[][], double) fill;
    Tuple!(double, double) stormEvent;
    Tuple!(double[], int) drainageSelect;
    int[] col;
    int[] row;
    int[] path;
    double[] z; 
    
    Grid grid;                      
    Cell[] cell;
    Sediment[] seds;
    Drainage[] drainageObj; 
    Stream stream;
    Discharge discharge;
    
    // read all input files and create objects
    seds = ReadSeds();                              // read sediment properties
    auto cellProps = ReadCells();                   // read cell elevations
    col = cellProps[0];
    row = cellProps[1];
    z = cellProps[2];
    grid = new Grid(col, row, z, seds.length);      // model geometry, including cell objects
    stream = new Stream(grid.d);                    // hydraulics
    discharge = new Discharge();                    // stochastic discharges
    drainageObj = ReadDrainage(grid);               // locations and other characteristics of drainages
    
    // initialize drainage random walk-based sediment compositions
    foreach (drainageItem; drainageObj){ 
        drainageSelect = drainageItem.Seed();
        drainageFrac ~= drainageSelect[0];          // note: drainageFrac and sedSelect are separate arrays not included in drainageObj
        sedSelect ~= drainageSelect[1];
        }
    
    dSed.length = seds.length;                      // initialize temporary array for handling deposition output
    dSed[] = 0.;
    
    // open discharge tabulation file
    auto dischargeFile = File("dischargeHistory.csv", "w"); 
    dischargeLine = "t" ~ "," ~ "Q" ~ "," ~ "dt";   
    dischargeFile.writeln(dischargeLine);
    
    // progress forward in time 
    while (t < discharge.endTime){          // outer loop = entire simulation, by discharge event
    
        // generate flow rate and time span (i.e., discharge event)
        stormEvent = discharge.Event();                     
        Q = stormEvent[0];
        dtEvent = stormEvent[1];            // dtEvent = time span associated with fixed discharge event
        dtEvent = FixTime(dtEvent, t, discharge.endTime);       // constrain dtEvent so that it does imply exceedance of endTime
        dischargeLine = to!string(t) ~ "," ~ to!string(Q) ~ "," ~ to!string(dtEvent);   
        dischargeFile.writeln(dischargeLine);           
    
        // loop through fan formation by drainage
        for (int k = 0; k < drainageObj.length; ++k){
    
            // update the drainage area sediment fraction composition
            drainageSelect = drainageObj[k].RandomWalk(sedSelect[k], drainageFrac[k]);  
            drainageFrac[k] = drainageSelect[0];
            sedSelect[k] = drainageSelect[1];
            
            // calculate path
            path = grid.Path(drainageObj[k].disCellIndex, drainageObj[k].edgeSrcIndex);         // return array of cell indices comprising fan rivulet
            Cell[] reach;
            foreach (indexNo; path){reach ~= grid.cell[indexNo];}
        
            // select debris flow or settling deposition regime
            Qeff = max(Q*drainageObj[k].fracQ, discharge.dryCutoffQ);
            if ((Qeff >= discharge.debrisFlowQ) && (uniform01() > 0.5)) {
                // debris flow
                J = stream.DebrisFlow(reach, Qeff, drainageFrac[k]);   // deposited thickness; J[numCells][numSeds]
                for (int i = 0; i < reach.length; ++i){
                    reach[i].Deposit(J[i]);}
                writeln("\tDebris flow occurred.");         
                }
            else{
                // deposition from stream flow
                tStep = 0.;             // elapsed time tracker during this fixed discharge event   
                while (tStep < dtEvent) {           // inner loop = inner time steps within discharge event (by change in sediment thickness)
                    // calculate deposition for this event
                    fill = stream.SettlingEvent(reach, Qeff, seds, drainageFrac[k]);            // note: placeholder to expand to multiple drainages
                    J = fill[0];            // deposition rates; J[numCells][numSeds]
                    dtStep = fill[1];       // time step associated with deposition rate (amount of time for specified accumulation)
                    dtStep = FixTime(dtStep, tStep, dtEvent);       // constrain dtStep so that tStep does not exceed dtEvent
                    // update sediment stacks in all cells in reach
                    for (int i = 0; i < reach.length; ++i){
                        for (int j = 0; j < seds.length; ++j){dSed[j] = J[i][j]*dtStep;}
                        if (sum(dSed)>0){reach[i].Deposit(dSed);}
                        }
                    tStep += dtStep;
                    }
                }

        }   // end of drainage-sweep loop
    
        t += dtEvent;
        writeln("Time = " ~ "\t" ~ to!string(t) ~ "\t" ~  "Q = " ~ "\t" ~ to!string(Q));
        
        }
        
    dischargeFile.close();                          // close the discharge tabulation file
    grid.WriteSedOuput(seds);                       // write sediments output file
    grid.WriteElevOuput(seds);                      // write cell land-surface elevations output file   
    writeln("Finished.");
    
    }