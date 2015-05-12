#ifndef Imagestats_h
#define Imagestats_h


#include <G3D/G3DAll.h>

class ImageStats
{
public:
	ImageStats(void);
	virtual ~ImageStats(void);
	
	typedef enum		{ FT_COLOR, FT_NORMAL, FT_TEXTURE, FT_DEPTH, FT_SIZE } feature;
	typedef enum		{ RED , GREEN, BLUE} channel;
	typedef	enum		{ MEAN, VAR} featureStat;

	void clear();
	void init(int width, int height, int sampleBudget, int samplesPerIteration);
	void addSample(int x, int y, Vector3& color,Vector3& normal,Vector3& texture,Vector3& depth);

	float getMean(int x, int y, int ft, int c){return m_dataMn[ft][c][x][y];}
	Color3 getMean(int x, int y, int ft){return Color3(m_dataMn[ft][RED][x][y],m_dataMn[ft][GREEN][x][y],m_dataMn[ft][BLUE][x][y]);}

	void setNewMean(int x, int y, float rMean, float gMean,float bMean);

	float getVariance(int x, int y, int ft, int c);
	Color3 getVariance(int x, int y, int ft);

	int width() {return m_w;}
	int height() {return m_h;}
	int getSamplePerPixel(int x, int y) {return m_samples[x][y];}
	int getSamplePerPixelToGo(int x, int y) {return m_samplesToGo[x][y];}
	float getSamplesPerPixelAvg();
	int getSampleBudget() { return m_sampleBudget;}
	int getSamplePerIteration() {return m_samplesPerIteration;}
	void substractFromSampleBudget(int samplesUsed) { m_sampleBudget -= samplesUsed;}
	void setAdaptiveSampling(int x, int y, int samplesToGo) {m_samplesToGo[x][y] = samplesToGo;}

	shared_ptr<Image3>& getImageMean(int ft, bool normalize = false);
	shared_ptr<Image3>& getImageVar(int ft, bool normalize = false);
	shared_ptr<Image3>& getImageMeanGammaCorrected(int ft, float gamma);
	
	
private:
	const static int COLOR_CHS = 3;

	float**				m_dataMn[FT_SIZE][COLOR_CHS]; //Data Mean for each channel
	float**				m_dataVr[FT_SIZE][COLOR_CHS]; //Data variance for each channel
	int**				m_samples; //number of samples per pixel
	int**				m_samplesToGo; //number of samples to render for adaptive sampling
	shared_ptr<Image3>	m_imageMean[FT_SIZE];
	shared_ptr<Image3>	m_imageVar[FT_SIZE];

	float				m_meanMin[FT_SIZE][COLOR_CHS];
	float				m_meanMax[FT_SIZE][COLOR_CHS];
	float				m_varMin[FT_SIZE][COLOR_CHS];
	float				m_varMax[FT_SIZE][COLOR_CHS];

	int		m_w;
	int		m_h;
	int		m_sampleBudget;
	int		m_samplesPerIteration;
	int		m_avgSamplesPP;
};

#endif