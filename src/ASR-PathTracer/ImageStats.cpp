#include "ImageStats.h"


ImageStats::ImageStats(void) :
	m_avgSamplesPP(1)
{
	for (int f = 0 ; f < FT_SIZE ; f++){
		for (int c = 0 ; c < COLOR_CHS ; c++){
			m_meanMin[f][c] = std::numeric_limits<float>::max();
			m_varMin[f][c] = std::numeric_limits<float>::max();
			m_meanMax[f][c] = std::numeric_limits<float>::min();
			m_varMax[f][c] = std::numeric_limits<float>::min();

			m_dataMn[f][c] = nullptr;
			m_dataVr[f][c] = nullptr;
		}
	}
	m_samples = nullptr;
	m_samplesToGo = nullptr;
}


ImageStats::~ImageStats(void)
{
	if (m_dataMn[0] != nullptr){
		for (int x = 0 ; x < m_w ; x++){
			delete m_samples[x];
			delete m_samplesToGo[x];
			for (int f = 0 ; f < FT_SIZE ; f++){
				delete m_dataMn[f][0][x]; delete m_dataMn[f][1][x]; delete m_dataMn[f][2][x];
				delete m_dataVr[f][0][x]; delete m_dataVr[f][1][x]; delete m_dataVr[f][2][x];
			}
		}
		delete m_samples;
		delete m_samplesToGo;
		for (int f = 0 ; f < FT_SIZE ; f++){
			delete m_dataMn[f][0]; delete m_dataMn[f][1]; delete m_dataMn[f][2];
			delete m_dataVr[f][0]; delete m_dataVr[f][1]; delete m_dataVr[f][2];
		}
	}
}

void ImageStats::init(int width, int height, int sampleBudget, int samplesPerIteration){

	for (int f = 0 ; f < FT_SIZE ; f++){
		if (!notNull(m_imageMean) || m_w != width || m_h != height)
			m_imageMean[f]	= Image3::createEmpty(width, height);
		if (!notNull(m_imageVar) || m_w != width || m_h != height)
			m_imageVar[f]	= Image3::createEmpty(width, height);
	}

	for (int f= 0 ; f < FT_SIZE ; f++){
		for (int c = 0 ; c < COLOR_CHS ; c++){ 
			if (m_dataMn[f][c] != nullptr) {
				for (int x = 0 ; x < m_w ; x++){
					delete m_dataMn[f][c][x];
					delete m_dataVr[f][c][x];
				}
				delete m_dataMn[f][c];
				delete m_dataVr[f][c];
			}
			m_dataMn[f][c] = new float*[width];
			m_dataVr[f][c] = new float*[width];
		}
	}

	m_samples		= new int*[width];
	m_samplesToGo	= new int*[width];
	for (int x = 0 ; x < width ; x++){
		for (int f = 0 ; f < FT_SIZE ; f++){
			m_dataMn[f][0][x] = new float[height];
			m_dataMn[f][1][x] = new float[height];
			m_dataMn[f][2][x] = new float[height];

			m_dataVr[f][0][x] = new float[height];
			m_dataVr[f][1][x] = new float[height];
			m_dataVr[f][2][x] = new float[height];
		}

		m_samples[x]		= new int[height];
		m_samplesToGo[x]	= new int[height];
		for (int y = 0 ; y < height ; y++){
			for (int f = 0 ; f < FT_SIZE ; f++){
				m_dataMn[f][0][x][y] = 0;
				m_dataMn[f][1][x][y] = 0;
				m_dataMn[f][2][x][y] = 0;

				m_dataVr[f][0][x][y] = 0;
				m_dataVr[f][1][x][y] = 0;
				m_dataVr[f][2][x][y] = 0;
			}

			m_samples[x][y] = 0;
			m_samplesToGo[x][y] = samplesPerIteration;
		}
	}
	m_w = width;
	m_h = height;
	m_sampleBudget = sampleBudget;
	m_samplesPerIteration = samplesPerIteration;
}

void ImageStats::clear()
{
	if (m_dataMn[0] != nullptr) {
		for (int x = 0 ; x < m_w ; x++)
			for (int y = 0 ; y < m_h ; y++){
				for (int f = 0 ; f < FT_SIZE ; f++){
					m_dataMn[f][0][x][y] = 0; m_dataMn[f][1][x][y] = 0; m_dataMn[f][2][x][y] = 0;
					m_dataVr[f][0][x][y] = 0; m_dataVr[f][1][x][y] = 0; m_dataVr[f][2][x][y] = 0;
				}
				m_samples[x][y] = 0;
			}
	}
}

//If there is just one value per image (like for depth) we use only the RED channel
void ImageStats::addSample(int x, int y, Vector3& color,Vector3& normal,Vector3& texture,Vector3& depth){
	float oldMean, newMean;
	m_samples[x][y]++;
	for (int c = 0 ; c < COLOR_CHS ; c++)
	{
		oldMean = m_dataMn[FT_COLOR][c][x][y];
		newMean = oldMean + (color[c] - oldMean)/m_samples[x][y];
		m_dataVr[FT_COLOR][c][x][y] += (color[c] - oldMean)*(color[c] - newMean); 
		m_dataMn[FT_COLOR][c][x][y] = newMean;

		oldMean = m_dataMn[FT_NORMAL][c][x][y];
		newMean = oldMean + (normal[c] - oldMean)/m_samples[x][y];
		m_dataVr[FT_NORMAL][c][x][y] += (normal[c] - oldMean)*(normal[c] - newMean); 
		m_dataMn[FT_NORMAL][c][x][y] = newMean;

		oldMean = m_dataMn[FT_TEXTURE][c][x][y];
		newMean = oldMean + (texture[c] - oldMean)/m_samples[x][y];
		m_dataVr[FT_TEXTURE][c][x][y] += (texture[c] - oldMean)*(texture[c] - newMean); 
		m_dataMn[FT_TEXTURE][c][x][y] = newMean;

		oldMean = m_dataMn[FT_DEPTH][c][x][y];
		newMean = oldMean + (depth[c] - oldMean)/m_samples[x][y];
		m_dataVr[FT_DEPTH][c][x][y] += (depth[c] - oldMean)*(depth[c] - newMean); 
		m_dataMn[FT_DEPTH][c][x][y] = newMean;
	}
	m_avgSamplesPP = 0;
}

float ImageStats::getVariance(int x, int y, int ft, int c){
	return (m_samples[x][y] < 2 ) ? 0.0 : m_dataVr[ft][c][x][y]/float(m_samples[x][y] - 1);
}
Color3 ImageStats::getVariance(int x, int y, int ft){
	return (m_samples[x][y] < 2 ) ? Color3(0.0) : Color3(m_dataVr[ft][RED  ][x][y]/float(m_samples[x][y] - 1),
														 m_dataVr[ft][GREEN][x][y]/float(m_samples[x][y] - 1),
														 m_dataVr[ft][BLUE ][x][y]/float(m_samples[x][y] - 1));
}

void ImageStats::setNewMean(int x, int y, float rMean, float gMean,float bMean){ 
	m_dataMn[FT_COLOR][RED][x][y] = rMean; 
	m_dataMn[FT_COLOR][GREEN][x][y] = gMean; 
	m_dataMn[FT_COLOR][BLUE][x][y] = bMean;
}

shared_ptr<Image3>& ImageStats::getImageMean(int ft, bool normalize){
	if (normalize){
		for (int c = 0 ; c < COLOR_CHS ; c++){
			m_meanMin[ft][c] = std::numeric_limits<float>::max();
			m_meanMax[ft][c] = std::numeric_limits<float>::min();
			for (int x = 0 ; x < m_w ; x++){
				for (int y = 0 ; y < m_h ; y++){
					if ( m_dataMn[ft][c][x][y] < m_meanMin[ft][c]) m_meanMin[ft][c] = m_dataMn[ft][c][x][y];
					if ( m_dataMn[ft][c][x][y] > m_meanMax[ft][c]) m_meanMax[ft][c] = m_dataMn[ft][c][x][y];
				}
			}
		}
		for (int x = 0 ; x < m_w ; x++)
			for (int y = 0 ; y < m_h ; y++)
				m_imageMean[ft]->fastSet(x,y,Color3((m_dataMn[ft][RED][x][y] - m_meanMin[ft][RED])/(m_meanMax[ft][RED] - m_meanMin[ft][RED]),
												(m_dataMn[ft][GREEN][x][y] - m_meanMin[ft][GREEN])/(m_meanMax[ft][GREEN] - m_meanMin[ft][GREEN]),
												(m_dataMn[ft][BLUE][x][y] - m_meanMin[ft][BLUE])/(m_meanMax[ft][BLUE] - m_meanMin[ft][BLUE])));
	}
	else{
		for (int x = 0 ; x < m_w ; x++)
			for (int y = 0 ; y < m_h ; y++)
				m_imageMean[ft]->fastSet(x,y,Color3(m_dataMn[ft][RED][x][y],m_dataMn[ft][GREEN][x][y],m_dataMn[ft][BLUE][x][y]));
	}
	return m_imageMean[ft];
}

shared_ptr<Image3>& ImageStats::getImageMeanGammaCorrected(int ft, float gamma){

	for (int x = 0 ; x < m_w ; x++)
		for (int y = 0 ; y < m_h ; y++)
			m_imageMean[ft]->fastSet(x,y,Color3(pow(m_dataMn[ft][RED][x][y],gamma),
												pow(m_dataMn[ft][GREEN][x][y],gamma),
												pow(m_dataMn[ft][BLUE][x][y],gamma)));
	return m_imageMean[ft];
}

shared_ptr<Image3>& ImageStats::getImageVar(int ft, bool normalize){
	if (normalize){
		for (int c = 0 ; c <  COLOR_CHS ; c++){
			m_varMin[ft][c] = std::numeric_limits<float>::max();
			m_varMax[ft][c] = std::numeric_limits<float>::min();
			for (int x = 0 ; x < m_w ; x++){
				for (int y = 0 ; y < m_h ; y++){
					if ( getVariance(x,y,ft,c) < m_varMin[ft][c])  m_varMin[ft][c] = getVariance(x,y,ft,c);
					if ( getVariance(x,y,ft,c) > m_varMax[ft][c])  m_varMax[ft][c] = getVariance(x,y,ft,c);
				}
			}
		}
		for (int x = 0 ; x < m_w ; x++)
			for (int y = 0 ; y < m_h ; y++)
				m_imageVar[ft]->fastSet(x,y,Color3((getVariance(x,y,ft,RED) - m_varMin[ft][RED])/(m_varMax[ft][RED] - m_varMin[ft][RED]),
											       (getVariance(x,y,ft,GREEN) - m_varMin[ft][GREEN])/(m_varMax[ft][GREEN] - m_varMin[ft][GREEN]),
											       (getVariance(x,y,ft,BLUE) - m_varMin[ft][BLUE])/(m_varMax[ft][BLUE] - m_varMin[ft][BLUE])));
	}
	else{
		for (int x = 0 ; x < m_w ; x++)
			for (int y = 0 ; y < m_h ; y++)
				m_imageVar[ft]->fastSet(x,y,Color3(getVariance(x,y,ft,RED),getVariance(x,y,ft,GREEN),getVariance(x,y,ft,BLUE)));
	}
	return m_imageVar[ft];
}

float ImageStats::getSamplesPerPixelAvg(){
	if (m_avgSamplesPP) return m_avgSamplesPP; //If already calculated just return them (whenever we add a sample we return this to 0 again)

	float count = 0;
	for(int x = 0 ; x < m_w ;x++)
		for(int y = 0 ; y < m_h ; y++)
			count+= m_samples[x][y];
	return count/float(m_w*m_h);
}
