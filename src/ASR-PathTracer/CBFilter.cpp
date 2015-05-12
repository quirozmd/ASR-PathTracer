#include "CBFilter.h"
#include <math.h> 

CBFilter::CBFilter(void) : 
	m_sigma_r(1.0),
	m_sigma_fn(0.8),
	m_sigma_ft(0.25),
	m_sigma_fd(0.6),
	m_sigmaPreSUREOpt(8),
	m_sumS(0)
{
	m_Samp = nullptr;
	for (int c = 0 ; c< COLOR_CHS ; c++){
		for (int s = 0 ; s < KERNEL_SIZES ; s++){
			m_F_ci[c][s]		= nullptr; 
			m_SUREerror[c][s]	= nullptr; 
			m_filterSUREerror[c][s]	= nullptr; 
		}
	}
	//Different sigma sizes for the filter bank
	if (KERNEL_SIZES == 7){
		kernelSigma[0] = 1;
		kernelSigma[1] = sqrt(2);
		kernelSigma[2] = 2;
		kernelSigma[3] = 2*sqrt(2);
		kernelSigma[4] = 4;
		kernelSigma[5] = 4*sqrt(2);
		kernelSigma[6] = 8;
	}
	else if (KERNEL_SIZES == 4){
		kernelSigma[0] = 1;
		kernelSigma[1] = 2;
		kernelSigma[2] = 4;
		kernelSigma[3] = 8;
	}

	m_sigma_f22[FT_COLOR]	= 2.0 * m_sigma_r * m_sigma_r;
	m_sigma_f22[FT_NORMAL]  = 2.0 * m_sigma_fn * m_sigma_fn;
	m_sigma_f22[FT_TEXTURE] = 2.0 * m_sigma_ft * m_sigma_ft;
	m_sigma_f22[FT_DEPTH]   = 2.0 * m_sigma_fd * m_sigma_fd;
}

CBFilter::~CBFilter(void)
{
	if ( m_Samp == nullptr) return;

	for (int x = 0 ; x < m_w ; x++){
		delete m_Samp[x];
		for (int c = 0 ; c< COLOR_CHS ; c++)
			for (int s = 0 ; s < KERNEL_SIZES ; s++){		
				delete m_F_ci[c][s][x];
				delete m_SUREerror[c][s][x];
				delete m_filterSUREerror[c][s][x];
			}
	}
	for (int c = 0 ; c< COLOR_CHS ; c++)
		for (int s = 0 ; s < KERNEL_SIZES ; s++){
			delete m_SUREerror[c][s];
			delete m_filterSUREerror[c][s];
			delete m_F_ci[c][s];
		}
}

void CBFilter::init(ImageStats* dataStats){

	m_dataStats = dataStats;
	m_w = dataStats->width();
	m_h = dataStats->height();

	m_optimizedImage = Image3::createEmpty(m_w, m_h);
	m_scaleSelectionMapImage = Image3::createEmpty(m_w, m_h);

	//Initialize data for each pixel in each channel
	for (int c = 0 ; c< COLOR_CHS ; c++){
		for (int s = 0 ; s < KERNEL_SIZES ; s++){
			m_SUREerror[c][s]		= new float*[m_w]; 
			m_filterSUREerror[c][s]	= new float*[m_w]; 
			m_F_ci[c][s]			= new float*[m_w]; 
		}
	}
	m_Samp = new float*[m_w];
	for (int x = 0 ; x < m_w ; x++){
		m_Samp[x] = new float[m_h];
		for (int c = 0 ; c< COLOR_CHS ; c++)
			for (int s = 0 ; s < KERNEL_SIZES ; s++){
				m_SUREerror[c][s][x]		= new float[m_h]; 
				m_filterSUREerror[c][s][x]	= new float[m_h]; 
				m_F_ci[c][s][x]				= new float[m_h]; 
			}
		for (int y = 0 ; y < m_h; y++){
			m_Samp[x][y] = 0.0;
			for (int c = 0 ; c< COLOR_CHS ; c++){
				for (int s = 0 ; s < KERNEL_SIZES ; s++){
					m_F_ci[c][s][x][y]				= 0.0; 
					m_SUREerror[c][s][x][y]			= 0.0; 
					m_filterSUREerror[c][s][x][y]	= 0.0; 
				}
			}
		}
	}

}

shared_ptr<Image3> CBFilter::applyFilter(int id){
	
	m_currentFilterSizeID = id;
	m_filteredImage[id]	= Image3::createEmpty(m_w, m_h);
	errorMap[id]		= Image3::createEmpty(m_w, m_h);
	filterErrorMap[id]	= Image3::createEmpty(m_w, m_h);
	GThread::runConcurrently2D(Point2int32(0, 0), Point2int32(m_w, m_h), this, &CBFilter::filterPixelWithError);

	//Set texture for later display
	char textureName[30];
	sprintf(textureName,"Error for sigma:%2.2f",kernelSigma[id]);
	errorMapText[id] = Texture::fromImage(textureName,errorMap[id]);
	sprintf(textureName,"Filtered sigma:%2.2f",kernelSigma[id]);
	m_filteredImageText[id] = Texture::fromImage(textureName,m_filteredImage[id]);

	return m_filteredImage[id];
}

shared_ptr<Image3> CBFilter::optimize(){

	GThread::runConcurrently2D(Point2int32(0, 0), Point2int32(m_w, m_h), this, &CBFilter::findMinErorAndOptimizePixel);

	m_scaleSelectionMapText = Texture::fromImage("Scale selection map",m_scaleSelectionMapImage);
	for (int s = 0 ; s < KERNEL_SIZES ; s++){
		char textureName[30];
		sprintf(textureName,"Filtered Error sigma:%2.2f",kernelSigma[s]);
		filterErrorMapText[s] = Texture::fromImage(textureName,filterErrorMap[s]);
	}

	return m_optimizedImage;
}

shared_ptr<Image3> CBFilter::adaptiveSampling(){
	m_samplingDensityImage = Image3::createEmpty(m_w, m_h);

	m_sumS = 0;
	for(int x = 0 ; x < m_w ;x++)
		for(int y = 0 ; y < m_h ; y++)
			m_sumS+= m_Samp[x][y];

	GThread::runConcurrently2D(Point2int32(0, 0), Point2int32(m_w, m_h), this, &CBFilter::adpativeSamplingToPixel);

	return m_samplingDensityImage;
}

void CBFilter::filterPixelWithError(int x, int y, int threadID){
	int X,Y;
	float w_ij[3], spacialTerm, featureTerm[FT_SIZE];
	float c_i[3]; //The color of the current pixel pix[x][y]
	float c_j[3]; //The color of the neighbour pixel pix[X][Y]
	float derivF; //The derivative of the filter dF(c_i)/dc_i
	float sigmaCi2[3]; //Squared Variance of the sample
	float m_weightSum[3];
	float m_F2_ci[3];
	float f_i_j[FT_SIZE][3],divisor;
	float sigma_i[FT_SIZE][3], sigma_j[FT_SIZE];
	float sigma_s22;
		
	sigma_s22 = 2.0 * kernelSigma[m_currentFilterSizeID] * kernelSigma[m_currentFilterSizeID];
	float size = ceil(2.0*kernelSigma[m_currentFilterSizeID]);

	for(int c = 0 ; c < COLOR_CHS ; c++){
		c_i[c] = m_dataStats->getMean(x,y,FT_COLOR,c);
		sigmaCi2[c] = m_dataStats->getVariance(x,y,FT_COLOR,c)/m_dataStats->getSamplesPerPixelAvg();

		m_F2_ci[c] = 0;
		m_F_ci[c][m_currentFilterSizeID][x][y] = 0;
		m_weightSum[c] = 0;
		
		for (int f= 0 ; f < FT_SIZE ; f++)
			sigma_i[f][c] = m_dataStats->getVariance(x,y,f,c);
	}
	
	//For each pixel wee need to find its neighbourhood defined by m_size
	for ( int nx = -size ; nx <= size ; nx++){
		for ( int ny = -size ; ny <= size ; ny++){
			X = x + nx; 
			Y = y + ny;
			if (X >= 0 && X < m_w && Y >= 0 && Y < m_h){
				spacialTerm = ( (X-x)*(X-x) + (Y-y)*(Y-y) ) / sigma_s22;
				for (int c = 0 ; c< COLOR_CHS ; c++){
					c_j[c] = m_dataStats->getMean(X,Y,FT_COLOR,c);
					for (int f = 0 ; f < FT_SIZE ; f++){
						f_i_j[f][c] = m_dataStats->getMean(X,Y,f,c) - m_dataStats->getMean(x,y,f,c);

						sigma_j[f] = m_dataStats->getVariance(X,Y,f,c);

						divisor = abs(sigma_i[f][c]+sigma_j[f]) + 0.0000001;
						featureTerm[f] =  (f_i_j[f][c]*f_i_j[f][c]) / (divisor*m_sigma_f22[f]);
						}
					w_ij[c] =	exp(-spacialTerm)*
							//exp(-featureTerm[FT_COLOR])*
							exp(-featureTerm[FT_DEPTH])*
							exp(-featureTerm[FT_TEXTURE])*
							exp(-featureTerm[FT_NORMAL]);
					m_weightSum[c] += w_ij[c];
				}
					
				for (int c = 0 ; c< COLOR_CHS ; c++){
					m_F_ci[c][m_currentFilterSizeID][x][y] += ((w_ij[0]+w_ij[1]+w_ij[2])/3.0)*c_j[c];
					m_F2_ci[c] += ((w_ij[0]+w_ij[1]+w_ij[2])/3.0)*c_j[c]*c_j[c];
				}
			}
		}
	}
	for (int c = 0 ; c< COLOR_CHS ; c++){
		m_F_ci[c][m_currentFilterSizeID][x][y] = m_F_ci[c][m_currentFilterSizeID][x][y]/((m_weightSum[0]+m_weightSum[1]+m_weightSum[2])/3.0);
		m_F2_ci[c] = m_F2_ci[c]/((m_weightSum[0]+m_weightSum[1]+m_weightSum[2])/3.0);

		//Sure error calculation before filtering
		derivF = (1.0/((m_weightSum[0]+m_weightSum[1]+m_weightSum[2])/3.0)) + ((m_F2_ci[c] - (m_F_ci[c][m_currentFilterSizeID][x][y]*m_F_ci[c][m_currentFilterSizeID][x][y]) )/(m_sigma_r*m_sigma_r));
		m_SUREerror[c][m_currentFilterSizeID][x][y] = ((m_F_ci[c][m_currentFilterSizeID][x][y] - c_i[c])*(m_F_ci[c][m_currentFilterSizeID][x][y] - c_i[c])) + (2.0*sigmaCi2[c]*derivF) - sigmaCi2[c];
	}
	m_filteredImage[m_currentFilterSizeID]->fastSet(x,y,Color3(m_F_ci[0][m_currentFilterSizeID][x][y],m_F_ci[1][m_currentFilterSizeID][x][y],m_F_ci[2][m_currentFilterSizeID][x][y]));
	Vector3 errorGrayScale(m_SUREerror[0][m_currentFilterSizeID][x][y],m_SUREerror[1][m_currentFilterSizeID][x][y],m_SUREerror[2][m_currentFilterSizeID][x][y]);
	errorMap[m_currentFilterSizeID]->fastSet(x,y,Color3(pow(errorGrayScale.squaredLength(),1.0/4.0)));
}

void CBFilter::findMinErorAndOptimizePixel(int x, int y, int threadID){
	////Filter error with another cross-bilateral filter to reduce variance
	int X,Y;
	float m_weightSum[3][KERNEL_SIZES];
	float sigma_s22,sigma_i[3][FT_SIZE],sigma_j[FT_SIZE];
	float w_ij[3], spacialTerm, rangeTerm, featureTerm[FT_SIZE],divisor;
	float c_i_j[3],f_i_j[FT_SIZE][3];

	for (int c = 0 ; c < COLOR_CHS ; c++){
		for (int f= 0 ; f < FT_SIZE ; f++)
			sigma_i[f][c] = m_dataStats->getVariance(x,y,f,c);
	}
	

	for (int s = 0 ; s < KERNEL_SIZES ; s++){

		for (int c = 0 ; c < COLOR_CHS ; c++){
			m_weightSum[c][s] = 0;
			m_filterSUREerror[c][s][x][y] = 0;
		}
		
		sigma_s22 = 2.0 * m_sigmaPreSUREOpt * m_sigmaPreSUREOpt;
		float size = ceil(2.0*m_sigmaPreSUREOpt);

		//For each pixel wee need to find its neighbourhood defined by m_size
		for ( int nx = -size ; nx <= size ; nx++){
			for ( int ny = -size ; ny <= size ; ny++){
				X = x + nx; 
				Y = y + ny;
				if (X >= 0 && X < m_w && Y >= 0 && Y < m_h){
					spacialTerm = ( (X-x)*(X-x) + (Y-y)*(Y-y) ) / sigma_s22;

					c_i_j[0] = m_SUREerror[0][s][X][Y] - m_SUREerror[0][s][x][y];
					c_i_j[1] = m_SUREerror[1][s][X][Y] - m_SUREerror[1][s][x][y];
					c_i_j[2] = m_SUREerror[2][s][X][Y] - m_SUREerror[2][s][x][y];

					rangeTerm = ( (c_i_j[0] * c_i_j[0]) + (c_i_j[1] * c_i_j[1]) + (c_i_j[2] * c_i_j[2]) ) / (m_sigma_f22[FT_COLOR]);

					for (int c = 0 ; c< COLOR_CHS ; c++){
						for (int f = 1 ; f < FT_SIZE ; f++){
							f_i_j[f][c] = m_dataStats->getMean(X,Y,f,c) - m_dataStats->getMean(x,y,f,c);

							sigma_j[f] = m_dataStats->getVariance(X,Y,f,c);

							divisor = abs(sigma_i[f][c]+sigma_j[f]) + 0.0000001;
							featureTerm[f] =  (f_i_j[f][c]*f_i_j[f][c]) / (divisor*m_sigma_f22[f]);
							}
						w_ij[c] =	exp(-spacialTerm)*
									//exp(-rangeTerm)*
									exp(-featureTerm[FT_DEPTH])*
									exp(-featureTerm[FT_TEXTURE])*
									exp(-featureTerm[FT_NORMAL]);
						m_weightSum[c][s] += w_ij[c];
					}		
					
					for (int c = 0 ; c< COLOR_CHS ; c++)
						m_filterSUREerror[c][s][x][y] += ((w_ij[0]+w_ij[1]+w_ij[2])/3.0)*m_SUREerror[c][s][X][Y];
				}
			}
		}
		for (int c = 0 ; c< COLOR_CHS ; c++)
			m_filterSUREerror[c][s][x][y] = m_filterSUREerror[c][s][x][y]/((m_weightSum[0][s]+m_weightSum[1][s]+m_weightSum[2][s])/3.0);

		Vector3 filtErrorGrayScale(m_filterSUREerror[0][s][x][y],m_filterSUREerror[1][s][x][y],m_filterSUREerror[2][s][x][y]);
		filterErrorMap[s]->fastSet(x,y,Color3(pow(filtErrorGrayScale.squaredLength(),1.0/4.0)));
	}
	

	int scaleSelection;
	float minError = std::numeric_limits<float>::max();
	for (int s = 0 ; s < KERNEL_SIZES ; s++){
		float error = (m_filterSUREerror[0][s][x][y]*m_filterSUREerror[0][s][x][y] + 
						m_filterSUREerror[1][s][x][y]*m_filterSUREerror[1][s][x][y] + 
						m_filterSUREerror[2][s][x][y]*m_filterSUREerror[2][s][x][y]);
		if (error < minError){
			minError = error;
			scaleSelection = s;
		}
	}
	m_dataStats->setNewMean(x,y,m_F_ci[RED][scaleSelection][x][y],
							    m_F_ci[GREEN][scaleSelection][x][y],
								m_F_ci[BLUE][scaleSelection][x][y]);

	m_scaleSelectionMapImage->fastSet(x,y,Color3((float(scaleSelection)/float(KERNEL_SIZES-1))));
	m_optimizedImage->fastSet(x,y,m_dataStats->getMean(x,y,FT_COLOR));

	//SURE(F(c_i))
	float SURE_F_i = sqrt((m_SUREerror[RED][scaleSelection][x][y]*m_SUREerror[RED][scaleSelection][x][y]) +
					      (m_SUREerror[GREEN][scaleSelection][x][y]*m_SUREerror[GREEN][scaleSelection][x][y]) +
				          (m_SUREerror[BLUE][scaleSelection][x][y]*m_SUREerror[BLUE][scaleSelection][x][y]));
	float sigma_i_2 = sqrt((m_dataStats->getVariance(x,y,FT_COLOR,0)*m_dataStats->getVariance(x,y,FT_COLOR,0)) +
				           (m_dataStats->getVariance(x,y,FT_COLOR,1)*m_dataStats->getVariance(x,y,FT_COLOR,1)) +
					       (m_dataStats->getVariance(x,y,FT_COLOR,2)*m_dataStats->getVariance(x,y,FT_COLOR,2))) ;
	float Luminance_i =	(0.2126f * m_F_ci[RED][scaleSelection][x][y]) + 
						(0.7152f * m_F_ci[GREEN][scaleSelection][x][y]) + 
						(0.0722  * m_F_ci[BLUE][scaleSelection][x][y]);
	//Set adaptive sampling function S(i)
	m_Samp[x][y] =  (SURE_F_i + sigma_i_2) / ( (Luminance_i*Luminance_i) + 0.001);
}

void CBFilter::adpativeSamplingToPixel(int x, int y, int threadID){
	if (m_sumS){
		float samplesToGo = m_w*m_h*m_Samp[x][y] / m_sumS;
		m_dataStats->setAdaptiveSampling(x,y,ceil(m_dataStats->getSamplePerIteration()*samplesToGo));
		m_samplingDensityImage->fastSet(x,y,Color3(samplesToGo/5.0));
	}
	else{
		m_dataStats->setAdaptiveSampling(x,y,m_dataStats->getSamplePerIteration());
		m_samplingDensityImage->fastSet(x,y,Color3(1.0));
	}
}
