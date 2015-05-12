#ifndef App_h
#define App_h

#include <G3D/G3DAll.h>
#include <GLG3D/GLG3D.h>
#include "ImageStats.h"
#include "CBFilter.h"
#include <fstream>
#include <iostream>
#include "World.h"

class App : public GApp {
private:
    
    int                 m_maxBounces;
	const float			BUMP_DISTANCE;
	float				m_image_scale;
	int                 m_sampleBudget;
	int					m_samplesPerIteration;
	int					m_sampleCount;
	Pointer<String>		m_sampleCountPtrString;
	String				m_sampleCountString;
	float				m_outputGammCorrection;
	int					m_imgWidth;
	int					m_imgHeight;
	World::Scene		m_worldScene;
	bool				m_readyToSample;
	bool				m_readyToFilter;
	bool				m_readyToOptimize;
	bool				m_debugMode;

	typedef	enum		{ MEAN, VAR} featureStat;
	typedef enum      	{ DISPLAY_RESULT,DISPLAY_TEXTURE } displayMode; //Display the result of path tracing 0 or just the texture 1(to easily place camera)

	std::ofstream		m_timeOutputFile;   //File to output times of each step of the algorithm
	Stopwatch			m_timer;		    //Timer to time each step.
	double				m_runSeconds;		//To sum up the whole running time withouth aiting time;
	void				tick();							//Starts the m_timer
	double				tock(const std::string &operation);			//Stops the m_timer and outputs in m_timeOutputFile the seconds for the last task, also add-up to m_runSeconds
	void				debugTime(const std::string &str);

	displayMode			m_displayMode;

    World*              m_world;

    Array<Random>       m_rng;

    /** Allocated by expose and render */
    shared_ptr<Texture> m_result;

	/** Array to store the name of each feature, mainly for GUI */
	String				m_featureNames[CBFilter::FT_SIZE]; 

	/** Textures to show each feature for the filter */
	shared_ptr<Texture> m_featureTextures[CBFilter::FT_SIZE][2];

	Plane				m_focalPlane;
	/** Image and statistics to compute and storage each feature of the scene:	0-Color (the result of the path Tracing algorithm), 
																				1-Normal (the normal of the fisrt hit object packed in a RGB),
																				2-Texture (Lambertian color or texture of the first hit surface, similar to ambient term), 
																				3-Depth (Depth of the first hit surface)*/
	ImageStats			m_featureData; 

	//Cross-bilateral filter, for some particular size, with its associated SURE-error
	CBFilter			m_CBFilter;

	GuiTextureBox*		m_textureBox[CBFilter::FT_SIZE][2];
	GuiTabPane*			m_featureTabPane;

	String				m_stillsPath;

    /** Position during the previous frame */
    CFrame              m_prevCFrame;

    /** Called from onInit() */
    void makeGUI();

    /** Trace a whole image. */
    void pathTraceImage();

	/** Trace a single path into the world, it sets each of the last 4 parameters as the result*/
    void tracePath(const Ray& ray, World* world, Random& rng , int px, int py);

	/** Sample one pixel of m_currentImage. Called on multiple threads. */
    void samplePixel(int x, int y, int threadID);

	/** Generates a random sample according to some defined distrubution return the distrubution's pdf */
	float generateSample(Vector3& normal, Random& rng, Vector3& wDirection);

	Color3 sampleLightArea( World* world, Random& rng,const Ray& ray, shared_ptr<Surfel> surfel);

    /** Show a full-screen message */
    void message(const String& msg) const;

	GuiTabPane* makeFeaturePane( GuiPane* p, featureStat stat, int controlWidth = 300);

	void updateAndSave(const std::string &str);
	


public:

    App(const GApp::Settings& settings , World::Scene scene, int spi, int tspp);

    virtual void onInit();

    /** Callback for the More samples button */
    void onAddSamples();

	/** Callback for the reset button */
    void onReset();

	/** Callback for the filter button */
    void onFilter();

	/** Callback for the SURE-Optimization button */
    void onOptimize();

    virtual void onGraphics(RenderDevice* rd, Array<shared_ptr<Surface> >& posed3D, Array<shared_ptr<Surface2D> >& posed2D);
    virtual void onCleanup();
};

#endif
