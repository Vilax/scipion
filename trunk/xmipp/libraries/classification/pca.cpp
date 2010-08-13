/***************************************************************************
 *
 * Authors:     Jorge Garcia de la Nava Ruiz (gdl@ac.uma.es)
 *              Carlos Oscar S. Sorzano (coss@cnb.csic.es)
 *              Alberto Pascual Montano (pascual@cnb.csic.es)
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/

#include <strstream>

#ifdef __sun
#include <ieeefp.h>
#endif

#include "pca.h"
#include "distance.h"

/**
* Calculate the eigenval/vecs
* Parameter: ts The vectors.
*/
void PCAAnalyzer::reset(ClassicTrainingVectors const &ts)
{
    std::vector<unsigned> dummy;
    for (unsigned i = 0;i < ts.size();i++)
        dummy.push_back(i);
    reset(ts, dummy);
}

/**
* Calculate the eigenval/vecs
* Parameter: ts The vectors.
* Parameter: idx The indexes of the vectors to use
*/
void PCAAnalyzer::reset(ClassicTrainingVectors const &ts, std::vector<unsigned> const & idx)
{
    std::vector<FeatureVector> a;
    int n = ts.dimension();
    a.resize(n);
    int verbosity = listener->getVerbosity();

    {
        if (verbosity)
            listener->OnReportOperation((std::string) "Normalizing....\n");
        if (verbosity == 1)
            listener->OnInitOperation(n);

        //Get the mean of the given cluster of vectors
        for (int k = 0;k < n;k++)
        {
            a[k].resize(n);
            Feature sum = 0.0;
            int l = 0;
            for (std::vector<unsigned>::const_iterator i = idx.begin();i != idx.end();i++)
            {
                if (finite(ts.itemAt(*i)[k]))
                {
                    sum += ts.itemAt(*i)[k];
                    l++;
                }
            }
            mean.push_back(sum / l);
            if (verbosity == 1)
                listener->OnProgress(k);
        }
        if (verbosity == 1)
            listener->OnProgress(n);

        if (verbosity == 1)
            listener->OnInitOperation(n);
        for (int i = 0;i < n;i++)
        {
            for (int j = 0;j <= i;j++)
            {
                Feature sum = 0.0;
                int l = 0;
                for (std::vector<unsigned>::const_iterator it = idx.begin();it != idx.end();it++)
                {
                    Feature d1 = ts.itemAt(*it)[i] - mean[i];
                    Feature d2 = ts.itemAt(*it)[j] - mean[j];
                    if (finite(d1) && finite(d2))
                    {
                        sum += d1 * d2;
                        l++;
                    }
                }
                if (l)
                    a[i][j] = a[j][i] = sum / l;
                else
                    a[i][j] = a[j][i] = 0;
            }
            if (verbosity == 1)
                listener->OnProgress(i);
        }
        if (verbosity == 1)
            listener->OnProgress(n);

        //  for(int i=0;i<n;i++)
        //   std::cout << a[i] << std::endl;
    }

    eigenval.resize(n);
    eigenvec.resize(n);
    set_Dimension(n);

    FeatureVector b;
    b.resize(n);
    FeatureVector z;
    z.resize(n);
    FeatureVector &d = eigenval;
    std::vector<FeatureVector> &v = eigenvec;

    for (int i = 0;i < n;i++)
    {
        v[i].resize(n);
        v[i][i] = 1.0;
        b[i] = d[i] = a[i][i];
    }

    int nrot = 0;

    if (verbosity)
        listener->OnReportOperation((std::string) "Diagonalizing matrix....\n");
    if (verbosity == 1)
        listener->OnInitOperation(50);

    //Jacobi method (it=iterationn number)
    for (int it = 1;it <= 50;it++)
    {
        if ((verbosity == 1) && (it == 1))
            listener->OnProgress(0);

        Feature tresh;
        Feature sm = 0.0;
        for (int ip = 0; ip < n - 1; ip++)
        {
            for (int iq = ip + 1; iq < n; iq++)
                sm += fabs(a[iq][ip]);
        }
        if (sm == 0.0)
        {//Done. Sort vectors
            for (int i = 0; i < n - 1; i++)
            {
                int k = i;
                Feature p = d[i];

                for (int j = i + 1; j < n; j++)
                    if (d[j] >= p)
                        p = d[k = j];

                if (k != i)
                {//Swap i<->k
                    d[k] = d[i];
                    d[i] = p;
                    FeatureVector t = v[i];
                    v[i] = v[k];
                    v[k] = t;
                }
            }
            if (verbosity == 1)
                listener->OnProgress(50);
            return;
        }

        if (it < 4)
            tresh = 0.2 * sm / (n * n);
        else
            tresh = 0;

        for (int ip = 0; ip < n - 1; ip++)
        {
            for (int iq = ip + 1; iq < n; iq++)
            {
                Feature g = 100.0 * fabs(a[iq][ip]);

                if (it > 4
                    && fabs(d[ip]) + g == fabs(d[ip])
                    && fabs(d[iq]) + g == fabs(d[iq]))
                    a[iq][ip] = 0.0;
                else if (fabs(a[iq][ip]) > tresh)
                {
                    Feature tau, t, s, c;
                    Feature h = d[iq] - d[ip];
                    if (fabs(h) + g == fabs(h))
                        t = a[iq][ip] / h;
                    else
                    {
                        Feature theta = 0.5 * h / a[iq][ip];
                        t = 1.0 / (fabs(theta) + sqrt(1.0 + theta * theta));
                        if (theta < 0.0)
                            t = -t;
                    }
                    c = 1.0 / sqrt(1 + t * t);
                    s = t * c;
                    tau = s / (1.0 + c);
                    h = t * a[iq][ip];
                    z[ip] -= h;
                    z[iq] += h;
                    d[ip] -= h;
                    d[iq] += h;
                    a[iq][ip] = 0.0;

#define rotate(a,i,j,k,l) \
    g = a[i][j]; \
    h = a[k][l]; \
    a[i][j] = g - s *(h + g*tau); \
    a[k][l] = h + s*(g - h*tau);

                    int j;
                    for (j = 0; j < ip; j++)
                    {
                        rotate(a, ip, j, iq, j)
                    }
                    for (j = ip + 1; j < iq; j++)
                    {
                        rotate(a, j, ip, iq, j)
                    }
                    for (j = iq + 1; j < n; j++)
                    {
                        rotate(a, j, ip, j, iq)
                    }
                    for (j = 0; j < n; j++)
                    {
                        rotate(v, ip, j, iq, j)
                    }

                    nrot += 1;
                }//if
            }//for iq
        }//for ip
        for (int ip = 0; ip < n; ip++)
        {
            b[ip] += z[ip];
            d[ip] = b[ip];
            z[ip] = 0.0;
        }

        if (verbosity == 1)
            listener->OnProgress(it - 1);

    }//for it

    if (verbosity == 1)
        listener->OnProgress(50);


    REPORT_ERROR(1, "too many Jacobi iterations");
}

/* Prepare for correlation ------------------------------------------------- */
void PCAAnalyzer::prepare_for_correlation()
{
    int nmax = D;
    int dmax = mean.size();

    // Initialize
    prod_ei_mean.resize(nmax);
    prod_ei_ei.resize(nmax);
    avg_ei.resize(nmax);

    // Compute products <ei,ei>, <ei,mean>
    for (int n = 0; n < nmax; n++)
    {
        prod_ei_ei[n] = prod_ei_mean[n] = avg_ei[n] = 0;
        for (int d = 0; d < dmax; d++)
        {
            prod_ei_mean[n] += eigenvec[n][d] * mean[d];
            prod_ei_ei[n] += eigenvec[n][d] * eigenvec[n][d];
            avg_ei[n] += eigenvec[n][d];
        }
        avg_ei[n] /= dmax;
    }

    // Compute product <mean,mean>
    prod_mean_mean = avg_mean = 0;
    for (int d = 0; d < dmax; d++)
    {
        prod_mean_mean += mean[d] * mean[d];
        avg_mean += mean[d];
    }
    avg_mean /= dmax;
}

/**Set identity matrix as eigenvector matrix*/
void PCAAnalyzer::setIdentity(int n)
{
    if (n < 0)
        n = 0;
    eigenval.resize(n);
    fill(eigenval.begin(), eigenval.end(), 1.0);
    eigenvec.resize(n);
    for (int i = 0;i < n;i++)
    {
        eigenvec[i].resize(n);
        fill(eigenvec[i].begin(), eigenvec[i].end(), 0.0);
        eigenvec[i][i] = 1.0;
    }
}

/* Components for variance ------------------------------------------------- */
int PCAAnalyzer::Dimension_for_variance(double th_var)
{
    int imax = eigenval.size();
    double sum = 0;
    th_var /= 100;
    for (int i = 0; i < imax; i++)
        sum += eigenval[i];

    double explained = 0;
    int i = 0;
    do
    {
        explained += eigenval[i++];
    }
    while (explained / sum < th_var);
    return i;
}

/* Project ----------------------------------------------------------------- */
void PCAAnalyzer::Project(FeatureVector &input, FeatureVector &output)
{
    if (input.size() != eigenvec[0].size())
        REPORT_ERROR(1, "PCA_project: vectors are not of the same size");

    int size = input.size();
    output.resize(D);
    for (int i = 0; i < D; i++)
    {
        output[i] = 0;
        // Comput the dot product between the input and the PCA vector[i]
        for (int j = 0; j < size; j++)
            output[i] += input[j] * eigenvec[i][j];
    }
}

/* Clear ------------------------------------------------------------------- */
void PCAAnalyzer::clear()
{
    set_Dimension(0);
    mean.clear();
    eigenvec.clear();
    eigenval.clear();
}

/* Show/read PCA ----------------------------------------------------------- */
std::ostream& operator << (std::ostream &out, const PCAAnalyzer &PC)
{
    out << "Relevant Dimension: " << PC.get_Dimension() << std::endl;
    out << "Mean vector: ";
    int size = PC.mean.size();
    out << "(" << size << ") ---> ";
    for (int j = 0; j < size; j++)
        out << PC.mean[j] << " ";
    out << std::endl;
    for (int i = 0; i < PC.get_Dimension(); i++)
    {
        out << PC.eigenval[i] << " (" << size << ") ---> ";
        for (int j = 0; j < size; j++)
            out << PC.eigenvec[i][j] << " ";
        out << std::endl;
    }
    return out;
}

std::istream& operator >> (std::istream &in, PCAAnalyzer &PC)
{
    PC.clear();
    int D;
    std::string read_line;
    getline(in, read_line);
    sscanf(read_line.c_str(), "Relevant Dimension: %d", &D);
    PC.set_Dimension(D);
    PC.eigenval.resize(D);
    PC.eigenvec.resize(D);

    int size;
    getline(in, read_line);
    sscanf(read_line.c_str(), "Mean vector: (%d) --->", &size);
    read_line.erase(0, read_line.find('>') + 1); // remove until --->
    PC.mean.resize(size);
    std::istrstream istr1(read_line.c_str());
    for (int j = 0; j < size; j++)
        istr1 >> PC.mean[j];

    for (int i = 0; i < D; i++)
    {
        getline(in, read_line);
        float f;
        sscanf(read_line.c_str(), "%f (%d) ---> ", &f, &size);
        PC.eigenval[i] = f;
        read_line.erase(0, read_line.find('>') + 1); // remove until --->
        PC.eigenvec[i].resize(size);
        std::istrstream istr2(read_line.c_str());
        for (int j = 0; j < size; j++)
            istr2 >> PC.eigenvec[i][j];
    }

    return in;
}

/* PCA set destructor ------------------------------------------------------ */
PCA_set::~PCA_set()
{
    int imax = PCA.size();
    for (int i = 0; i < imax; i++)
        delete PCA[i];
}

/* Create empty PCA -------------------------------------------------------- */
int PCA_set::create_empty_PCA(int n)
{
    int retval = PCA.size();
    PCA.resize(retval + n);
    for (int i = 0; i < n; i++)
        PCA[retval+i] = new PCAAnalyzer;
    return retval;
}

/* Show/Read PCAset -------------------------------------------------------- */
std::ostream& operator << (std::ostream &out, const PCA_set &PS)
{
    int imax = PS.PCA.size();
    out << "Number of PCAs: " << imax << std::endl;
    for (int i = 0; i < imax; i++)
        out << *(PS.PCA[i]);
    return out;
}

std::istream& operator >> (std::istream &in, PCA_set &PS)
{
    int imax;
    std::string read_line;
    getline(in, read_line);
    sscanf(read_line.c_str(), "Number of PCAs: %d\n", &imax);
    PS.PCA.resize(imax);
    for (int i = 0; i < imax; i++)
    {
        PS.PCA[i] = new PCAAnalyzer;
        in >> *(PS.PCA[i]);
    }
    return in;
}

/* Running PCA constructor ------------------------------------------------- */
Running_PCA::Running_PCA(int _J, int _d)
{
    J = _J;
    d = _d;
    n = 0;
    sum_all_samples.initZeros(d);
    current_sample_mean.initZeros(d);
    sum_proj.initZeros(d);
    sum_proj2.initZeros(d);
    eigenvectors.initZeros(d, J);
}

/* Update with new sample -------------------------------------------------- */
void Running_PCA::new_sample(const Matrix1D<double> &sample)
{
    n++;

    // Re-estimate sample mean
    sum_all_samples += sample;
    current_sample_mean = sum_all_samples / n;

    // Estimate eigenvectors
    Matrix1D<double> un = sample;
    un.row = false;
    for (int j = 0; j < J; j++)
    {
        if (n <= j + 1)
        {
            // If there are not enough samples to estimate this eigenvector, then
            // skip it
            double norm = un.module();
            if (norm > XMIPP_EQUAL_ACCURACY)
                un /= norm;
            eigenvectors.setCol(j, un);
        }
        else
        {
            // If there are enough samples
            // Substract the sample mean to have a zero-mean vector
            if (j == 0)
                un -= current_sample_mean;

            // Compute the scale of this vector as the dot product
            // between un and the current eigenvector estimate
            double scale = 0;
            for (int i = 0; i < d; i++)
                scale += un(i) * eigenvectors(i, j);

            // Re-estimate the eigenvector
            double norm2 = 0;
            double w1 = (double)(n - 1.0) / n;
            double w2 = (1.0 - w1) * scale;
            for (int i = 0; i < d; i++)
            {
                eigenvectors(i, j) =
                    w1 * eigenvectors(i, j) +
                    w2 * un(i);
                norm2 += eigenvectors(i, j) *
                         eigenvectors(i, j);
            }

            // Renormalize
            if (norm2 > XMIPP_EQUAL_ACCURACY)
            {
                double norm = sqrt(norm2);
                for (int i = 0; i < d; i++)
                    eigenvectors(i, j) /= norm;
            }

            // Project un onto the space spanned by this eigenvector
            double project = 0;
            for (int i = 0; i < d; i++)
                project += un(i) * eigenvectors(i, j);
            for (int i = 0; i < d; i++)
                un(i) -= project * eigenvectors(i, j);

            // Update the variance of this vector
            sum_proj(j) += project;
            sum_proj2(j) += project * project;
        }
    }
}

/* Project a sample vector on the PCA space -------------------------------- */
void Running_PCA::project(const Matrix1D<double> &input,
                          Matrix1D<double> &output) const
{
    output.initZeros(J);
    for (int j = 0; j < J; j++)
        for (int i = 0; i < d; i++)
            output(j) += (input(i) - current_sample_mean(i)) * eigenvectors(i, j);
}

/* Subtract average ------------------------------------------------------- */
void PCAMahalanobisAnalyzer::subtractAvg()
{
    int N=v.size();
    if (N==0)
        return;
    // Compute average
    MultidimArray<double> avg;
    typeCast(v[0],avg);
    for (int n=1; n<N; n++)
    {
        MultidimArray<float> &aux=v[n];
        FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(avg)
        DIRECT_A1D_ELEM(avg,i)+=DIRECT_A1D_ELEM(aux,i);
    }
    avg/=N;

    // Subtract average and compute stddev
    MultidimArray<float> avgF;
    typeCast(avg,avgF);
    for (int n=0; n<N; n++)
        v[n]-=avgF;
}

/* Standardize variables -------------------------------------------------- */
void PCAMahalanobisAnalyzer::standardarizeVariables()
{
    int N=v.size();
    if (N==0)
        return;
    // Compute average
    MultidimArray<double> avg;
    typeCast(v[0],avg);
    for (int n=1; n<N; n++)
    {
        MultidimArray<float> &aux=v[n];
        FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(avg)
        DIRECT_A1D_ELEM(avg,i)+=DIRECT_A1D_ELEM(aux,i);
    }
    avg/=N;

    // Subtract average and compute stddev
    MultidimArray<float> avgF;
    typeCast(avg,avgF);
    MultidimArray<double> stddev;
    stddev.initZeros(avg);
    for (int n=0; n<N; n++)
    {
        v[n]-=avgF;
        MultidimArray<float> &aux=v[n];
        FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(aux)
        {
            float f=DIRECT_A1D_ELEM(aux,i);
            DIRECT_A1D_ELEM(stddev,i)+=f*f;
        }
    }
    stddev/=N-1;
    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(stddev)
    DIRECT_A1D_ELEM(stddev,i)=sqrt(DIRECT_A1D_ELEM(stddev,i));

    // Divide by stddev
    MultidimArray<float> istddevF;
    istddevF.resize(stddev);
    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(stddev)
    DIRECT_A1D_ELEM(istddevF,i)=(float)(1.0/DIRECT_A1D_ELEM(stddev,i));
    for (int n=0; n<N; n++)
    {
        MultidimArray<float> &aux=v[n];
        FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(aux)
        DIRECT_A1D_ELEM(aux,i)*=DIRECT_A1D_ELEM(istddevF,i);
    }
}

/* Add vector ------------------------------------------------------------- */
void PCAMahalanobisAnalyzer::projectOnPCABasis(Matrix2D<double> &CtY)
{
    int N=v.size();
    int NPCA=PCAbasis.size();
    CtY.initZeros(NPCA,N);
    for (int ii=0; ii<N; ii++)
    {
        const MultidimArray<float> &Iii=v[ii];
        for (int jj=0; jj<NPCA; jj++)
        {
            const MultidimArray<double> &Ijj=PCAbasis[jj];
            CtY(jj,ii)=0;
            FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(Iii)
            MAT_ELEM(CtY,jj,ii)+=
                DIRECT_A1D_ELEM(Iii,i)*DIRECT_A1D_ELEM(Ijj,i);
        }
    }
}

void PCAMahalanobisAnalyzer::learnPCABasis(int NPCA, int Niter)
{
    // Take the first vectors for the PCA basis
    MultidimArray<double> vPCA;
    NPCA=XMIPP_MIN(NPCA,v.size());
    for (int n=0; n<NPCA; n++)
    {
        typeCast(v[n],vPCA);
        PCAbasis.push_back(vPCA);
    }

    int N=v.size();
    for (int n=0; n<Niter; n++)
    {
        // E-step ..........................................................
        // Compute C^t*C
        Matrix2D<double> CtC(NPCA,NPCA);
        for (int ii=0; ii<NPCA; ii++)
        {
            const MultidimArray<double> &Iii=PCAbasis[ii];
            for (int jj=ii; jj<NPCA; jj++)
            {
                const MultidimArray<double> &Ijj=PCAbasis[jj];
                CtC(ii,jj)=0;
                FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(Ijj)
                MAT_ELEM(CtC,ii,jj)+=
                    DIRECT_A1D_ELEM(Iii,i)*DIRECT_A1D_ELEM(Ijj,i);
                if (ii!=jj)
                    CtC(jj,ii)=CtC(ii,jj);
            }
        }

        // Compute C^t*Y
        Matrix2D<double> CtY, X;
        projectOnPCABasis(CtY);
        X=CtC.inv()*CtY;

        // M-step ..........................................................
        Matrix2D<double> XtXXtinv=X.transpose()*(X*X.transpose()).inv();
        for (int ii=0; ii<NPCA; ii++)
        {
            MultidimArray<double> &Ipca=PCAbasis[ii];
            Ipca.initZeros();
            for (int jj=0; jj<N; jj++)
            {
                const MultidimArray<float> &I=v[jj];
                double val=XtXXtinv(jj,ii);
                FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(I)
                DIRECT_A1D_ELEM(Ipca,i)+=DIRECT_A1D_ELEM(I,i)*val;
            }
        }
    }

    // Normalize output vectors
    for (int ii=0; ii<NPCA; ii++)
    {
        double norm=sqrt(PCAbasis[ii].sum2());
        PCAbasis[ii]/=norm;
    }
}

/* Compute Statistics ----------------------------------------------------- */
void PCAMahalanobisAnalyzer::computeStatistics(MultidimArray<double> & avg,
        MultidimArray<double> & stddev)
{
    int N=v.size();
    if (N==0)
    {
        avg.clear();
        stddev.clear();
        return;
    }
    typeCast(v[0],avg);
    MultidimArray<double> aux;
    for (int n=1; n<N; n++)
    {
        typeCast(v[n],aux);
        avg+=aux;
    }
    avg/=N;
    stddev.initZeros(avg);
    for (int n=0; n<N; n++)
    {
        typeCast(v[n],aux);
        aux-=avg;
        aux*=aux;
        stddev+=aux;
    }
    double iN_1=1.0/(N-1);
    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(stddev)
    DIRECT_A1D_ELEM(stddev,i)=sqrt(DIRECT_A1D_ELEM(stddev,i)*iN_1);
}

/* Evaluate score --------------------------------------------------------- */
//#define DEBUG
void PCAMahalanobisAnalyzer::evaluateZScore(int NPCA, int Niter)
{
    int N=v.size();
    if (N==1)
    {
        Zscore.initZeros(N);
        idx=Zscore.indexSort();
        return;
    }
    
#ifdef DEBUG
    std::cout << "Input vectors\n";
    for (int n=0; n<N; n++)
    {
        FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(v[n])
        std::cout << DIRECT_A1D_ELEM(v[n],i) << " ";
        std::cout << std::endl;
    }
#endif
    subtractAvg();
#ifdef DEBUG

    std::cout << "\n\nInput vectors after normalization\n";
    for (int n=0; n<N; n++)
    {
        FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(v[n])
        std::cout << DIRECT_A1D_ELEM(v[n],i) << " ";
        std::cout << std::endl;
    }
#endif

    learnPCABasis(NPCA, Niter);
#ifdef DEBUG

    std::cout << "\n\nPCA basis\n";
    for (int n=0; n<NPCA; n++)
    {
        FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(PCAbasis[n])
        std::cout << DIRECT_A1D_ELEM(PCAbasis[n],i) << " ";
        std::cout << std::endl;
    }
#endif

    Matrix2D<double> proj;
    projectOnPCABasis(proj);

#ifdef DEBUG

    std::cout << "\n\nPCA projected values\n" << proj << std::endl;
#endif

    // Estimate covariance matrix
    Matrix2D<double> cov(NPCA,NPCA);
    for (int ii=0; ii<N; ii++)
    {
        for (int i=0; i<NPCA; i++)
            for (int j=0; j<NPCA; j++)
                cov(i,j)+=proj(i,ii)*proj(j,ii);
    }
    cov/=N;
#ifdef DEBUG

    std::cout << "\n\nCovariance matrix\n" << cov << std::endl;
#endif

    Matrix2D<double> covinv=cov.inv();
    Zscore.initZeros(N);
    for (int ii=0; ii<N; ii++)
    {
        Zscore(ii)=0;
        for (int i=0; i<NPCA; i++)
        {
            double xi=MAT_ELEM(proj,i,ii);
            for (int j=0; j<NPCA; j++)
                DIRECT_A1D_ELEM(Zscore,ii)+=xi*MAT_ELEM(proj,j,ii)*MAT_ELEM(covinv,i,j);
        }
        Zscore(ii)=sqrt(ABS(Zscore(ii)));
    }
#ifdef DEBUG
    std::cout << "\n\nZscores\n";
    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(Zscore)
    std::cout << DIRECT_A1D_ELEM(Zscore,i) << " ";
    std::cout << std::endl;
#endif

    idx=Zscore.indexSort();
}
#undef DEBUG
