#ifndef Cbfilter_h
#define Cbfilter_h


#include "ImageStats.h"

class CBFilter
{
public:
	CBFilter(void);
	~CBFilter(void);
	void init(ImageStats* dataStats);
	shared_ptr<Image3> applyFilter(int id);
	shared_ptr<Image3> adaptiveSampling();
	shared_ptr<Image3> optimize();
	shared_ptr<Texture> getErrorTexture(int id) { return errorMapText[id]; }
	shared_ptr<Texture> getFilterImgTexture(int id) { return m_filteredImageText[id]; }
	shared_ptr<Texture> getFilterErrorTexture(int id) { return filterErrorMapText[id]; }
	shared_ptr<Texture> getscaleSelectionMapTexture() { return m_scaleSelectionMapText; }
	typedef enum		{ FT_COLOR, FT_NORMAL, FT_TEXTURE, FT_DEPTH, FT_SIZE } feature;
	typedef enum		{ RED , GREEN, BLUE} channel;
	int					getFilterBankSize(){return KERNEL_SIZES;}


private:

	/** Compute filtered image and SURE error. Called on multiple threads. */
    void filterPixelWithError(int x, int y, int threadID);
	/** Find min error between all filters in filter bank and optimize the filter applied. Called on multiple threads. */
	void findMinErorAndOptimizePixel(int x, int y, int threadID);
	/** Compute the new sampling rate for each pixel and set it in the dataStats. Called on multiple threads. */
	void adpativeSamplingToPixel(int x, int y, int threadID);

	const static int COLOR_CHS = 3;
	const static int KERNEL_SIZES = 4;
	float kernelSigma[KERNEL_SIZES];
	int m_currentFilterSizeID;

	float**	m_F_ci[COLOR_CHS][KERNEL_SIZES];	// F(ci) the filter
	float**	m_SUREerror[COLOR_CHS][KERNEL_SIZES]; // Stein's Unbiased Risk Estimator error
	float**	m_filterSUREerror[COLOR_CHS][KERNEL_SIZES]; // Stein's Unbiased Risk Estimator error
	float** m_Samp; //Sampling rate for each pixel S(i) 
	float   m_sumS; //The sum of the sampling function S(i) of all pixels

	int m_w;
	int m_h;
	const float m_sigmaPreSUREOpt; //Global scale sigma to reduce variance
	const float m_sigma_r;	//Sigma parameter for color
	const float m_sigma_fn; //Sigma parameter for normal
	const float m_sigma_ft; //Sigma parameter for texture
	const float m_sigma_fd; //Sigma parameter for depth
	float		m_sigma_f22[FT_SIZE];

	ImageStats* m_dataStats;

	shared_ptr<Image3> m_optimizedImage;
	shared_ptr<Image3> m_scaleSelectionMapImage;
	shared_ptr<Texture> m_scaleSelectionMapText;
	shared_ptr<Image3> m_samplingDensityImage;

	shared_ptr<Image3> m_filteredImage[KERNEL_SIZES];
	shared_ptr<Image3> errorMap[KERNEL_SIZES];
	shared_ptr<Image3> filterErrorMap[KERNEL_SIZES];
	shared_ptr<Texture> m_filteredImageText[KERNEL_SIZES];
	shared_ptr<Texture> errorMapText[KERNEL_SIZES];
	shared_ptr<Texture> filterErrorMapText[KERNEL_SIZES];
};

#endif

