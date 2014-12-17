/*==========================================================
 * hmmviterbi_topK_cond_table_c.cpp
 *
 *
 * The calling syntax is:
 *
 *		[condVitBound,condVitLogP,currentState,logP,condVitBound2,condVitLogP2] = hmmviterbi_topK_cond_table_cpp(seq,tr,e,k);
 *
 * This is a MEX-file for MATLAB.
 *
 *========================================================*/




#include <stdio.h>
#include "mex.h"
#include <vector>
#include "math.h"
#include "Eigen/Dense"
#include <algorithm>
using namespace Eigen;
using namespace std;
double SMALLNUM = 1e-10;
double BIGNUM = 1e100;


//void main(){
//    mxArray *plhs[10000];
//    const mxArray *prhs[4] = {{},{},{},{}};
//    mexFunction( 6, plhs,
//                4, prhs);
//}
void printMatrix(const MatrixXd &x);
pair<MatrixXd,MatrixXd> TopKInMatrix(const MatrixXd &mat,const unsigned char k);



/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    double *seqtemp;
    vector<unsigned char> seq;              /* input scalar */
    double *tr;               /* 1xN input matrix */
    double *e;
    unsigned char k;
    size_t ncols;                   /* size of matrix */
    double *outMatrix;              /* output matrix */
    
    /* check for proper number of arguments */
    if(nrhs!=4) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs","Two inputs required.");
    }
    /* make sure the first input argument is scalar */
    if( !mxIsDouble(prhs[0])) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notScalar","seq must be a int vector.");
    }
    if( !mxIsDouble(prhs[1])) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notScalar","tr must be a double matrix.");
    }
    if( !mxIsDouble(prhs[2])) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notScalar","e must be a double matrix.");
    }
    if( !mxIsDouble(prhs[3])||
        mxGetNumberOfElements(prhs[3])!=1) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notScalar","k must be a scalar.");
    }
    
    
    
    /* get the value of the scalar input  */
    //be careful
    seqtemp = mxGetPr(prhs[0]);
    
    size_t seq_size = mxGetNumberOfElements(prhs[0]);
    tr = mxGetPr(prhs[1]);
    size_t tr_size = mxGetNumberOfElements(prhs[1]);
    e = mxGetPr(prhs[2]);
    size_t e_size = mxGetNumberOfElements(prhs[2]);
    k = (unsigned char)mxGetScalar(prhs[3]);
  
    
    
    /* get dimensions of the input matrix */
    int seq_length = mxGetN(prhs[0]);
    int num_states = mxGetM(prhs[1]);
    int num_emissions = mxGetN(prhs[2]);
    
    
    for (int i=0; i<seq_length; i++) {
        seq.push_back(seqtemp[i]);
    }
    
    
    MatrixXd log_tr(num_states,num_states);
    MatrixXd log_e(num_states,num_emissions);

    for (int j=0; j<num_states; j++) {
        for (int i=0; i<num_states; i++) {
            log_tr(i,j) = logf(tr[i+num_states*j]);
        }
    }
    //printMatrix(log_tr);
    for (int j=0; j<num_emissions; j++){
        for (int i=0; i<num_states; i++){
            log_e(i,j) = logf(e[i+num_states*j]);
        }
    }
    //printMatrix(log_e);

    //  ptr: numstate L <2 k>
    vector<vector<MatrixXd> > pTR(num_states,vector<MatrixXd>(seq_length));
    for (int i=0; i<num_states; i++) {
        for (int j=0; j<seq_length; j++) {
            pTR[i][j] = MatrixXd::Zero(2,k);
        }
    }
    
    // vf: L <numstate k>
    vector<MatrixXd>  vf(seq_length);
    for (int i=0; i<seq_length; i++) {
        vf[i] = MatrixXd::Zero(num_states,k);
    }
    MatrixXd v(num_states,k);
    v.setConstant(-BIGNUM);
    v(0,0) = 0;
    MatrixXd v_old = v;
    
    //loop through

    
    for (int count = 0; count<seq_length; count++) {
        for (int state = 0; state<num_states; state++) {
            MatrixXd val = v_old + log_tr.col(state).replicate(1,k);
            //mexPrintf("%d,%d\n",count,state);
            pair<MatrixXd,MatrixXd> resPair = TopKInMatrix(val,k);
            MatrixXd &bestPTRK = resPair.first;
            MatrixXd &bestValK = resPair.second;

            pTR[state][count] = bestPTRK;
            MatrixXd temp(1,k);
            temp.setConstant(log_e(state,seq[count]-1));
//                            if(count == 0){
//                                //printMatrix(temp);
//                                mexPrintf("%d\n\n",seq[count]);
//                            }
            v.row(state) =  temp + bestValK;
        }
        vf[count] = v;
        v_old = v;
    }

    MatrixXd finalState(2,k);
    MatrixXd AlllogP(1,k);
    pair<MatrixXd,MatrixXd> resPair = TopKInMatrix(v,k);
    finalState = resPair.first;
    AlllogP = resPair.second;
    //printMatrix(finalState);
    //printMatrix(AlllogP);
    
    
    //backtrace
    // L 2 k
    vector<MatrixXd>  currentState(seq_length);
    for (int i=0; i<seq_length-1; i++) {
        currentState[i] = MatrixXd::Zero(2,k);
    }
    currentState[seq_length-1] = finalState;
    
    
    for (int count = seq_length-2; count>=0; count--) {
        for (int kk = 0; kk<k; kk++) {
            size_t thisState_r = int(currentState[count+1](0,kk));
            size_t thisK_c = int(currentState[count+1](1,kk));
            currentState[count].col(kk) = pTR[thisState_r][count+1].col(thisK_c);
        }
    }
    MatrixXd AllStateBound(1,k);
    AllStateBound.setZero();
    for (int kk = 0; kk<k; kk++) {
        for (int count = 0; count<seq_length; count++) {
            AllStateBound(0,kk) +=  double(currentState[count](0,kk)>SMALLNUM);
        }
    }
//    printMatrix(AllStateBound);
//    printMatrix(AlllogP);
//    

    
    //all state bound
    
    // currentState L 2 k
    for (int i=0; i<seq_length; i++) {
        currentState[i] = MatrixXd::Zero(2,k);
    }
    
    
    //////////////backward
    
    //  Bptr: numstate L <2 k>
    vector<vector<MatrixXd> > BpTR(num_states,vector<MatrixXd>(seq_length));
    for (int i=0; i<num_states; i++) {
        for (int j=0; j<seq_length; j++) {
            BpTR[i][j] = MatrixXd::Zero(2,k);
        }
    }
    
    v.setConstant(-BIGNUM);
    v.col(0) = log_e.col(seq[seq_length-1]-1);
    
    // vb: L <numstate k>
    vector<MatrixXd>  vb(seq_length);
    for (int i=0; i<seq_length-1; i++) {
        vb[i] = MatrixXd::Zero(num_states,k);
    }
    vb[seq_length-1] = v;
    v_old = v;
    
    for (int count = seq_length-2; count>=0; count--) {
        for (int state = 0; state<num_states; state++) {
            MatrixXd val = v_old.colwise() + log_tr.row(state).transpose();
            pair<MatrixXd,MatrixXd> resPair = TopKInMatrix(val,k);
            MatrixXd &bestPTRK = resPair.first;
            MatrixXd &bestValK = resPair.second;
            BpTR[state][count] = bestPTRK;
            MatrixXd temp(num_states,k);
            v.row(state) = temp.setConstant(log_e(state,seq[count]-1)) + bestValK;
        }
        vb[count] = v;
        v_old = v;
    }
    

    
    // condition on 1
    // numstate L 1 k
    vector<vector<MatrixXd> > nBoundtable(num_states,vector<MatrixXd>(seq_length));
    for (int i=0; i<num_states; i++) {
        for (int j=0; j<seq_length; j++) {
            nBoundtable[i][j] = MatrixXd::Zero(1,k);
        }
    }
    vector<vector<MatrixXd> > logPtable(num_states,vector<MatrixXd>(seq_length));
    for (int i=0; i<num_states; i++) {
        for (int j=0; j<seq_length; j++) {
            logPtable[i][j] = MatrixXd::Zero(1,k);
        }
    }
    
  
    
    
    for (int count = 0; count<seq_length; count++) {
        for (int state = 0; state<num_states; state++) {
            MatrixXd temp(2,k);
            for (int kk=0; kk<k; kk++) {
                temp(0,kk) = double(state);
                temp(1,kk) = double(kk);
            }
            currentState[count] = temp;
            for (int fcount = count-1; fcount>=0; fcount--) {
                for (int kk=0; kk<k; kk++) {
                    size_t thisState_r = int(currentState[fcount+1](0,kk));
                    size_t thisK_c = int(currentState[fcount+1](1,kk));
                    currentState[fcount].col(kk) = pTR[thisState_r][fcount+1].col(thisK_c);
                }

            }
            vector<MatrixXd> FcurrentState(currentState);
            MatrixXd FlogP(1,k);
            FlogP = vf[count].row(state);
            
            //backward
//           currentState[seq_length-1] = temp;
            for (int bcount = count+1; bcount<seq_length; bcount++) {
                for (int kk=0; kk<k; kk++) {
                    size_t thisState_r = int(currentState[bcount-1](0,kk));
                    size_t thisK_c = int(currentState[bcount-1](1,kk));
                    currentState[bcount].col(kk) = BpTR[thisState_r][bcount-1].col(thisK_c);
                }
                
            }
            vector<MatrixXd> BcurrentState(currentState);
            MatrixXd BlogP(1,k);
            BlogP = vb[count].row(state);
            //printMatrix(FlogP);
            //printMatrix(BlogP);
            //
            MatrixXd temp2(k,k);
            
            for (int k2=0; k2<k; k2++) {
                for (int k1=0; k1<k; k1++) {
                    temp2(k1,k2) = FlogP(0,k1)+BlogP(0,k2);
                }
            }
            
            //printMatrix(temp2);
            pair<MatrixXd,MatrixXd> resPair = TopKInMatrix(temp2,k);
            MatrixXd &bestPTRK = resPair.first;
            MatrixXd &bestValK = resPair.second;
            for (int kk = 0; kk<k; kk++) {
                for (int fcount = 0;fcount<count; fcount++) {
                    nBoundtable[state][count](0,kk) += double(FcurrentState[fcount](0,bestPTRK(0,kk))>SMALLNUM);
                }
                for (int bcount = count;bcount<seq_length; bcount++) {
                    nBoundtable[state][count](0,kk) += double(BcurrentState[bcount](0,bestPTRK(1,kk))>SMALLNUM);
                }
                logPtable[state][count](0,kk) = bestValK(0,kk) - log_e(state,seq[count]-1);
            }
        }
    }
    
//    printMatrix(nBoundtable[0][0]);
//    printMatrix(logPtable[0][0]);
    
    //condition on 2 table
    // L 1 k
    vector<MatrixXd>  nBoundtable2(seq_length);
    nBoundtable2[0] = nBoundtable[0][0];
    for (int i=1; i<seq_length; i++) {
        nBoundtable2[i] = MatrixXd::Zero(1,k);
    }
    
    vector<MatrixXd>  logPtable2(seq_length);
    logPtable2[0] = logPtable[0][0];
    for (int i=1; i<seq_length; i++) {
        logPtable2[i] = MatrixXd::Zero(1,k);
    }

    for (int count = 1; count<seq_length; count++) {
        MatrixXd temp(2,k);
        for (int kk=0; kk<k; kk++) {
            temp(0,kk) = 0;
            temp(1,kk) = double(kk);
        }
        currentState[count-1] = temp;
        for (int fcount = count-2; fcount>=0; fcount--) {
            for (int kk=0; kk<k; kk++) {
                size_t thisState_r = int(currentState[fcount+1](0,kk));
                size_t thisK_c = int(currentState[fcount+1](1,kk));
                currentState[fcount].col(kk) = pTR[thisState_r][fcount+1].col(thisK_c);
            }
        }
        vector<MatrixXd> FcurrentState(currentState);
        MatrixXd FlogP(1,k);
        FlogP = vf[count-1].row(0);
        
        //backward
        currentState[count] = temp;
        for (int bcount = count+1; bcount<seq_length; bcount++) {
            for (int kk=0; kk<k; kk++) {
                size_t thisState_r = int(currentState[bcount-1](0,kk));
                size_t thisK_c = int(currentState[bcount-1](1,kk));
                currentState[bcount].col(kk) = BpTR[thisState_r][bcount-1].col(thisK_c);
            }
            
        }
        vector<MatrixXd> BcurrentState(currentState);
        MatrixXd BlogP(1,k);
        BlogP = vb[count].row(0);
  
        MatrixXd temp2(k,k);
        
        for (int k2=0; k2<k; k2++) {
            for (int k1=0; k1<k; k1++) {
                temp2(k1,k2) = FlogP(0,k1)+BlogP(0,k2);
            }
        }
        
        //printMatrix(temp2);
        pair<MatrixXd,MatrixXd> resPair = TopKInMatrix(temp2,k);
        MatrixXd &bestPTRK = resPair.first;
        MatrixXd &bestValK = resPair.second;
        for (int kk = 0; kk<k; kk++) {
            for (int fcount = 0;fcount<count; fcount++) {
                nBoundtable2[count](0,kk) += double(FcurrentState[fcount](0,bestPTRK(0,kk))>SMALLNUM);
            }
            for (int bcount = count;bcount<seq_length; bcount++) {
                nBoundtable2[count](0,kk) += double(BcurrentState[bcount](0,bestPTRK(1,kk))>SMALLNUM);
            }
            logPtable2[count](0,kk) = bestValK(0,kk) - log_e(0,seq[count]-1);
        }
    }
    
//        printMatrix(nBoundtable2[0]);
//        printMatrix(logPtable2[0]);

    
    
    //output
    
    
    /* create the output matrix */
    const mwSize dims[3] = {k,num_states,seq_length};
    plhs[0] = mxCreateNumericArray(3,dims,mxDOUBLE_CLASS,mxREAL);
    plhs[1] = mxCreateNumericArray(3,dims,mxDOUBLE_CLASS,mxREAL);
    plhs[2] = mxCreateDoubleMatrix(1,k,mxREAL);
    plhs[3] = mxCreateDoubleMatrix(k,1,mxREAL);
    plhs[4] = mxCreateDoubleMatrix(seq_length,k,mxREAL);
    plhs[5] = mxCreateDoubleMatrix(seq_length,k,mxREAL);


    /* get a pointer to the real data in the output matrix */
    double *p1 = mxGetPr(plhs[0]);
    double *p2 = mxGetPr(plhs[1]);
    double *p3 = mxGetPr(plhs[2]);
    double *p4 = mxGetPr(plhs[3]);
    double *p5 = mxGetPr(plhs[4]);
    double *p6 = mxGetPr(plhs[5]);
    for (int kk=0; kk<k; kk++) {
        for (int i=0; i<num_states; i++) {
            for (int j=0; j<seq_length; j++) {
                p1[kk+k*i+k*num_states*j] = nBoundtable[i][j](0,kk);
                p2[kk+k*i+k*num_states*j] = logPtable[i][j](0,kk);
            }
        }
    }
    for (int kk=0; kk<k; kk++) {
        p3[kk] = AllStateBound(0,kk);
        p4[kk] = AlllogP(0,kk);
    }
    for (int i=0; i<seq_length; i++) {
        for (int kk=0; kk<k; kk++) {
            p5[i+seq_length*kk] = nBoundtable2[i](0,kk);
            p6[i+seq_length*kk] = logPtable2[i](0,kk);
        }
    }
   
    
// never do this:http://www.mathworks.com/matlabcentral/newsreader/view_thread/333650
//    mxFree(tr);
//    mxFree(e);
//    mxFree(seq);
    
    



}

// override the original log
double logf(double num){
    if(abs(num)<SMALLNUM){
        return -BIGNUM;
    }
    else{
        return log(num);
    }
}

void printMatrix(const MatrixXd &x){
    for (int i=0; i<x.rows();i++ ) {
        for (int j=0; j<x.cols(); j++) {
            if(x(i,j)<-1e10){
               mexPrintf("-inf\t");
               continue;
            }
            mexPrintf("%lf\t",x(i,j));
        }
        mexPrintf("\n");
    }
}

typedef pair<double,int> IdPair;
bool pair_greater ( const IdPair& l, const IdPair& r){
    return l.first > r.first;
}

pair<MatrixXd,MatrixXd> TopKInMatrix(const MatrixXd &mat,const unsigned char k){
    size_t m = mat.rows();
    size_t n = mat.cols();
    vector<IdPair> vec;
    for (int j=0; j<n; j++) {
        for (int i=0; i<m; i++) {
            vec.push_back(make_pair(mat(i,j),i+m*j));
        }
    }
    nth_element(vec.begin(), vec.begin() + k, vec.end(), pair_greater);
    MatrixXd bestPTRK(2,k);
    MatrixXd bestValK(1,k);
    for (int i=0; i<k; i++) {
        bestValK(0,i) = vec[i].first;
        bestPTRK(0,i) = double(int(vec[i].second)%m);
        bestPTRK(1,i) = double(int(vec[i].second)/m);
    }
    return make_pair(bestPTRK,bestValK);
}









