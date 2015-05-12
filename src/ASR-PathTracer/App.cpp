#include "App.h"
#include "World.h"

G3D_START_AT_MAIN();

int main(int argc, char** argv) {
    GApp::Settings settings;
    settings.window.caption     = "Adaptive Sampled Reconstructed Path Tracing";


	World::Scene scene = World::SIBENIK ;//or World::CRYTEK_SPONZA or World::CORNELL_BOX
	int samplesPerIteration = 16;
	int totalSamplesPerPixel = 1;
	int resolution = 512;
	if (argc > 3)
	{
		scene = World::Scene(atoi(argv[1]));
		resolution = atoi(argv[2]);
		samplesPerIteration = atoi(argv[3]);
		totalSamplesPerPixel = atoi(argv[4]);
	}

	settings.window.width       = resolution; 
    settings.window.height      = resolution;
	return App(settings,scene,samplesPerIteration,totalSamplesPerPixel).run();
}

App::App(const GApp::Settings& settings, World::Scene scene, int spi, int tspp) : 
		GApp(settings),
		m_image_scale(0.25),
		m_maxBounces(50), 
		m_sampleBudget(tspp),
		m_samplesPerIteration(spi),
		m_stillsPath("stills\\"),
		m_worldScene(scene),
		m_world(NULL),
		m_runSeconds(0),
		BUMP_DISTANCE(0.001f),
		m_sampleCount(0),
		m_displayMode(DISPLAY_RESULT),// or DISPLAY_TEXTURE
		m_debugMode(true)
		{
	m_outputGammCorrection = 1.0/2.2;

	m_featureNames[CBFilter::FT_COLOR] = String("Color");
	m_featureNames[CBFilter::FT_NORMAL] = String("Normal");
	m_featureNames[CBFilter::FT_TEXTURE] = String("Texture");
	m_featureNames[CBFilter::FT_DEPTH] = String("Depth");
	catchCommonExceptions = false;

	m_imgWidth  = int(window()->width()*m_image_scale );
	m_imgHeight = int(window()->height()*m_image_scale);

	m_readyToSample = true;
	m_readyToFilter = false;
	m_readyToOptimize = false;
}


void App::onInit() {
	tick();

    message("Loading...");
	renderDevice->setSwapBuffersAutomatically(true);

	m_world = new World(m_worldScene);

    // Create one random number generator per thread
    m_rng.resize(GThread::numCores());
    for (int i = 0; i < m_rng.size(); ++i) {
        //m_rng[i].reset(0xF018A4D2 ^ i, false);
		//m_rng[i].reset(uint32(System::time()) ^ i);
    }
	Array<Plane,10> clipPlanes;
	float a,b,c,d;

    showRenderingStats = false;
    createDeveloperHUD();
    developerWindow->setVisible(false);
    developerWindow->cameraControlWindow->setVisible(false);
    m_debugCamera->filmSettings().setAntialiasingEnabled(false);
	m_debugCamera->filmSettings().setAntialiasingHighQuality(false);
	m_debugCamera->filmSettings().setBloomStrength(0.0f);
	m_debugCamera->setFrame(m_world->frame);
	m_debugCamera->setFieldOfView(m_world->FOVDeg*pi()/180.0,FOVDirection::VERTICAL);
	m_debugCamera->getClipPlanes(Rect2D(Vector2(window()->width(),window()->height())),clipPlanes);

	clipPlanes[0].getEquation(a,b,c,d);
	m_focalPlane = Plane::fromEquation(-a,-b,-c,-(d-(m_world->focalDist)));	
	
    makeGUI();
    // Force re-render on first frame
    m_prevCFrame = CFrame(Matrix3::zero());
	float time = tock("Loading scene");

	//Initialize StatsData class 
	m_featureData.init(m_imgWidth,m_imgHeight,m_sampleBudget,m_samplesPerIteration);
	//Initialize filter with the features data pointers;
	m_CBFilter.init(&m_featureData);
}


void App::makeGUI() {
    shared_ptr<GuiWindow> window = GuiWindow::create("Controls", debugWindow->theme(), Rect2D::xywh(0,0,0,0), GuiTheme::TOOL_WINDOW_STYLE);
    GuiPane* pane = window->pane();
    
	//pane->addButton("Restart", this, &App::onReset);
	//if (m_debugMode){
		pane->addButton("Filter", this, &App::onFilter);
		pane->addButton("SURE-Optimization", this, &App::onOptimize);
		pane->addButton("More samples", this, &App::onAddSamples);
	//}
	m_sampleCountPtrString = Pointer<String>(&m_sampleCountString);
	pane->addTextBox("Samples: ",m_sampleCountPtrString,G3D::GuiTextBox::DELAYED_UPDATE,G3D::GuiTheme::NO_BACKGROUND_UNLESS_FOCUSED_TEXT_BOX_STYLE );
    
    //pane->addNumberBox("Samples ppx", &m_samplesPerPixel, "", GuiTheme::LINEAR_SLIDER, 1, 5000, 1);
    //pane->addNumberBox("Max bounces", &m_maxBounces, "", GuiTheme::LINEAR_SLIDER, 1, 30, 1);

	window->pack();
    window->setVisible(true);
    addWidget(window);

	debugPane->moveBy(200, 10);
    debugPane->setCaption(GuiText("\t\t Mean\t\t\t\t\t Variance", GFont::fromFile(System::findDataFile("arial.fnt")), 16));
	
	debugPane->setVisible(true);
	debugPane->beginRow();
	GuiTabPane* tabPane = dynamic_cast<GuiTabPane*>(makeFeaturePane(debugPane,MEAN));
	tabPane->pack();

	GuiTabPane* tabPane2 = dynamic_cast<GuiTabPane*>(makeFeaturePane(debugPane,VAR));
	tabPane2->pack();
	tabPane2->moveRightOf(tabPane,50);
	debugPane->endRow();

	debugPane->pack();   

	char numStr[20];
	sprintf(numStr,"%d",m_sampleCount);
	m_sampleCountPtrString.setValue(String(numStr));	

	
}


void App::onGraphics(RenderDevice* rd, Array<shared_ptr<Surface> >& surface3D, Array<shared_ptr<Surface2D> >& surface2D) {
        // Update the preview image only while moving
	if (m_debugMode){
		if (! m_prevCFrame.fuzzyEq(m_debugCamera->frame())) {
			onAddSamples();
			m_prevCFrame = m_debugCamera->frame();
		}
	}
	else{
		if (m_featureData.getSamplesPerPixelAvg() < m_sampleBudget){
			//onAddSamples();
			if (m_readyToSample)
				onAddSamples();
			else{
				if (m_readyToFilter)
					onFilter();
				else{
					if (m_readyToOptimize)
						onOptimize();
				}
			}
		}
		else{
			if (m_readyToFilter){
				onFilter();
				onOptimize();
			}
		}
	}


    if (m_result) {
        rd->push2D(); {
			if (m_displayMode == DISPLAY_RESULT)
				Draw::rect2D(rd->viewport(), rd, Color3::white(), m_result);
			else
				Draw::rect2D(rd->viewport(), rd, Color3::white(), m_featureTextures[CBFilter::FT_TEXTURE][MEAN]);
        } rd->pop2D();
    }

    Surface2D::sortAndRender(rd, surface2D);

}

void App::onReset() {
	for (int j = 0 ; j < CBFilter::FT_SIZE ; j++) 
		m_featureData.clear();

	m_sampleCount = 0;

	char numStr[20];
	sprintf(numStr,"%d",m_sampleCount);
	m_sampleCountPtrString.setValue(String(numStr));

	//onRender();
}

void App::onAddSamples() {
	//message("Sampling..."); 
	App::window()->setCaption("(Sampling...) Adaptive Sampled Reconstructed Path Tracing");
	shared_ptr<Texture> samplingDensity = Texture::fromImage("Sampling density",m_CBFilter.adaptiveSampling());
	//if (m_debugMode) 
	//	show(samplingDensity);
	pathTraceImage();
	m_readyToSample = false;
	m_readyToFilter = true;
	updateAndSave("");
}

void App::onFilter() {
    //message("Filtering..."); 
	App::window()->setCaption("(Filtering...) Adaptive Sampled Reconstructed Path Tracing");
	double totalFilterTime = 0;
	char numStr[50];
	shared_ptr<Image3> filteredFixSize;
	for (int s = 0 ; s < m_CBFilter.getFilterBankSize() ; s++){
		tick();
		filteredFixSize = m_CBFilter.applyFilter(s);
		if (m_debugMode) {
			show(m_CBFilter.getFilterImgTexture(s));
			show(m_CBFilter.getErrorTexture(s));
		}

		sprintf(numStr,"Filter id %d:",s);
		totalFilterTime += tock(numStr);
	}
	shared_ptr<Texture> filteredText = Texture::fromImage("Filtered img",filteredFixSize);
	m_film->exposeAndRender(renderDevice, m_debugCamera->filmSettings(), filteredText, m_result);

	sprintf(numStr,"Total filter time \t%4.2f",totalFilterTime);
	debugTime(numStr);

	m_readyToFilter = false;
	m_readyToOptimize = true;
}

void App::onOptimize() {
	tick();

    //message("SURE Optimization..."); 
	App::window()->setCaption("(SURE Optimization...) Adaptive Sampled Reconstructed Path Tracing");
	shared_ptr<Texture> filteredText = Texture::fromImage("Filtered img",m_CBFilter.optimize());
	m_film->exposeAndRender(renderDevice, m_debugCamera->filmSettings(), filteredText, m_result);
	if (m_debugMode){
		for (int s = 0 ; s < m_CBFilter.getFilterBankSize() ; s++)
			show(m_CBFilter.getFilterErrorTexture(s));
		show(m_CBFilter.getscaleSelectionMapTexture());
	}
	
	tock("Sure Optimization");
	m_readyToOptimize = false;
	m_readyToSample = true;
	updateAndSave("_Opt");
}

void App::updateAndSave(const std::string &str){	
	String mainFileName = m_stillsPath+m_world->sceneName;
	
	char widthNum[10],heightNum[10],samplesNum[10];
	float sampleNumber = m_featureData.getSamplesPerPixelAvg();
	sprintf(widthNum,"%d",m_imgWidth);
	sprintf(heightNum,"%d",m_imgHeight);
	sprintf(samplesNum,"%4.2f",sampleNumber);
	mainFileName.append("_" + String(widthNum) + "x" + String(heightNum) + "_" + samplesNum);
	
	m_featureData.getImageMeanGammaCorrected(CBFilter::FT_COLOR,m_outputGammCorrection)->save(mainFileName + str.c_str() + ".png");


	//if (sampleNumber > 1){
	//	FileSystem::createDirectory(mainFileName);
	//	for (int f = 0 ; f < CBFilter::FT_SIZE ; f++){
	//		m_featureData.getImageMean(f,f == CBFilter::FT_DEPTH)->save(mainFileName + "\\" + m_world->sceneName + "_" + m_featureNames[f] + "_MEAN" + ".png");
	//		m_featureData.getImageVar(f,false)->save(mainFileName + "\\" + m_world->sceneName + "_" + m_featureNames[f] + "_VAR" + ".png");
	//		if (f == CBFilter::FT_COLOR)
	//			m_featureData.getImageMean(f)->save(mainFileName + ".png");
	//	}
	//}
}

void App::pathTraceImage() {
	tick();    

	for (int i = 0; i < m_rng.size(); ++i)
		m_rng[i].reset(uint32(System::time()) ^ i);
	
	GThread::runConcurrently2D(Point2int32(0, 0), Point2int32(m_imgWidth, m_imgHeight), this, &App::samplePixel);
	
	float sampleAvg = m_featureData.getSamplesPerPixelAvg();
	m_featureData.substractFromSampleBudget(sampleAvg);

	
	char numStr[20];
	sprintf(numStr,"%2.2f",sampleAvg);
	m_sampleCountPtrString.setValue(String(numStr));

	

	// Post-processing
	shared_ptr<Texture> src = Texture::fromImage("Source", m_featureData.getImageMean(CBFilter::FT_COLOR));
	if (m_result) {
		m_result->resize(m_imgWidth, m_imgHeight);
	}
	
	m_film->exposeAndRender(renderDevice, m_debugCamera->filmSettings(), src, m_result);

	for (int j = 0; j < CBFilter::FT_SIZE; j++) {
		m_featureTextures[j][MEAN] = Texture::fromImage(m_featureNames[j]+String(" mean"),m_featureData.getImageMean(j,j == CBFilter::FT_DEPTH));
		m_featureTextures[j][VAR] = Texture::fromImage(m_featureNames[j]+String(" variance"),m_featureData.getImageVar(j,false));
		for (int k = 0 ; k < 2 ; k++)
		{
			m_textureBox[j][k]->setTexture(m_featureTextures[j][k]);
			m_textureBox[j][k]->zoomToFit();
		}	
	}
	
	tock("Path tracing");
}

void App::samplePixel(int x, int y, int threadID) {
	float randX,randY,randZ;
    Random& rng = m_rng[threadID];
	for (int i = 0; i < m_featureData.getSamplePerPixelToGo(x,y) ; ++i){
		rng.sphere(randX,randY,randZ);
		randX = (randX/2.0) + 0.5; 
		randY = (randY/2.0) + 0.5;
		Ray initialRay = m_debugCamera->worldRay(x + randX, y + randY ,Rect2D(Vector2(m_imgWidth,m_imgHeight)));
		Point3 pointInFocalPlane = initialRay.intersection(m_focalPlane);
		
		rng.sphere(randX,randY,randZ);
		Point3 newCamPos = initialRay.origin() + (Vector3(randX,randY,0.0)*m_world->blurRadius);
		Vector3 newDir = (pointInFocalPlane - newCamPos).unit();
		tracePath(Ray(newCamPos,newDir), m_world, rng,x,y);
	}
}

void App::tracePath(const Ray& ray, World* world, Random& rng , int px, int py) {
	int bounce = 0; //Current path length
	Ray r = ray;
    Color3 cl = Color3(0.0f,0.0f,0.0f); //acumulated color
	Color3 cf = Color3(1.0f,1.0f,1.0f); //acumulated reflectance

	Vector3 dataToAdd[CBFilter::FT_SIZE];

	while (1){
		float dist = std::numeric_limits<float>::max();
		const shared_ptr<Surfel>& surfel = world->intersect(r, dist);

		if (!notNull(surfel)){
			cl += m_world->ambient;
			break ;//If ray path hits nothing, return accumulated color plus the sky color (ambient term)
		}
		
		//Properties of the intersection
		Vector3 x  = surfel->position; 
		Vector3 n  = surfel->shadingNormal; // or shadingNormal
		Vector3 nl = n.dot(r.direction()) < 0 ? n : n*-1; //Properly oriented surface normal
		nl = nl.unit();
		Color3	f  = surfel->reflectivity(rng); //surfel->finiteScatteringDensity(w_i, -ray.direction())
		Color3  emi = surfel->emittedRadiance(-ray.direction());

		//Collect aditional data only for the first bounce
		if (bounce == 0){
			dataToAdd[CBFilter::FT_TEXTURE] = Vector3(f.r,f.g,f.b);
			dataToAdd[CBFilter::FT_NORMAL] = (nl + Vector3(1,1,1))*0.5;
			dataToAdd[CBFilter::FT_DEPTH] = x;
		}
		
		//Radiance from sampling Light (for Explict PT)
		Color3 LES = Color3(0.0,0.0,0.0);
		//if (emi.max() == 0.0)
		LES = sampleLightArea(world,rng,r,surfel);
		if (bounce > 0) emi = Color3(0.0);

		cl += cf*(emi+LES);

		// max reflectivity for Russian Roulette
		float p = f.max();
		bounce++;
		if (bounce > 5){  // Dont use russian roullette until 5th bounce
			if (rng.uniform() < p) 
				f = f*(1.0/p); 
			else 
				break ;
		}

		if (bounce >= m_maxBounces || m_displayMode == DISPLAY_TEXTURE)
			break ;

		Vector3 sampledDir;
		float pdf = generateSample(nl,rng,sampledDir);
		Color3 brdf = f*pdf;
		

		cf *= brdf;
		r = Ray::fromOriginAndDirection(x + sampledDir*BUMP_DISTANCE, sampledDir);
		debugAssert(r.direction().isFinite());
	}
	cl = Color3(clamp(cl.r,0.0f,1.0f),clamp(cl.g,0.0f,1.0f),clamp(cl.b,0.0f,1.0f));
	dataToAdd[CBFilter::FT_COLOR] = Vector3(cl.r,cl.g,cl.b);

	m_featureData.addSample(px,py,dataToAdd[CBFilter::FT_COLOR],dataToAdd[CBFilter::FT_NORMAL],dataToAdd[CBFilter::FT_TEXTURE],dataToAdd[CBFilter::FT_DEPTH]);
}

float App::generateSample(Vector3& normal,Random& rng, Vector3& wDirection){
	// x,y,z coordinates of the sample in the sphere
	float r1 = 2.0*pif()*rng.uniform(); 
    float r2 = sqrt(rng.uniform());

	float sampleX, sampleY, sampleZ;
	sampleX = cos(r1)*sqrt(1.0 - r2); 
    sampleY = sin(r1)*sqrt(1.0 - r2); 
    sampleZ = r2;
	Vector3 w = normal;
    // (w.y > 1 - eps) is equivalent to (w.x^2 + w.z^2) = ||w \cross (0,1,0)||^2 < 0.1
	Vector3 u = w.cross((abs(w.y) > 0.9) ? Vector3(1,0,0) : Vector3(0,1,0));
	u = u.unit();
    Vector3 v = w.cross(u);
	wDirection = Vector3(sampleX*u.x + sampleY*v.x + sampleZ*w.x,
						 sampleX*u.y + sampleY*v.y + sampleZ*w.y,
						 sampleX*u.z + sampleY*v.z + sampleZ*w.z);
	wDirection = wDirection.unit();
	//For convinience this pdf is actually cos(theta)/PI*pdf (to just multiply it by the diffuse color of the object)
	return 1.0; // real pdf: 1 / 4PI
}

Color3 App::sampleLightArea( World* world,Random& rng,const Ray& ray, shared_ptr<Surfel> surfel){
	Color3 e = Color3(0.0f,0.0f,0.0f);

	Sphere sp; //For each sphere light
	for (int L = 0 ; L < world->lightSurfacesIndex.size() ; L++ ) {
		world->getLightSpherei(L,sp);
		
		//Create coordinate system for sampling sw, su, sv
		Vector3 sw = sp.center - surfel->position; // vector in the direction of the light center
        // (sw.y > 1 - eps) is equivalent to (sw.x^2 + sw.z^2) = ||sw \cross (0,1,0)||^2 < 0.1
		Vector3 su = sw.cross((abs(sw.y) > 0.9) ? Vector3(1,0,0) : Vector3(0,1,0));
		su = su.unit();
        Vector3 sv = sw.cross(su);

        // Determine max angle
        float cos_a_max = sqrt(1.0 - sp.radius*sp.radius/sw.dot(sw));

        // Calculate sample direction (to light) based acording to equation from Realistic Rendering
		float eps1 = rng.uniform();
        float eps2 = rng.uniform();
        float cos_a = 1.0 - eps1 + (eps1*cos_a_max);
        float sin_a = sqrt(1.0 - cos_a*cos_a);
        float phi = 2.0*pif()*eps2;
        float xx = cos(phi)*sin_a;
        float yy = sin(phi)*sin_a;
        float zz = cos_a;
        Vector3 ld = Vector3(xx*su.x + yy*sv.x + zz*sw.x,
							 xx*su.y + yy*sv.y + zz*sw.y,
							 xx*su.z + yy*sv.z + zz*sw.z);
		ld = ld.unit();
        // Check for oclussion with shadow ray
		float dist = std::numeric_limits<float>::max();
		Vector3 nl = surfel->shadingNormal.dot(ray.direction()) < 0 ? surfel->shadingNormal : surfel->shadingNormal*-1;
		const shared_ptr<Surfel>& lightSurfel = world->intersect(Ray(surfel->position+nl*BUMP_DISTANCE,ld), dist);
		if ( notNull(lightSurfel) && lightSurfel->emittedRadiance(ld).max() > 0.0) {
			float omega = 2.0*pif()*(1.0 - cos_a_max);
			e += (surfel->reflectivity(rng)*(m_world->getEmissiveFromLighti(L)*ld.dot(nl)*omega))*(1.0/pif());
		}
		
	}
    return e;
}

void App::onCleanup() {
    delete m_world;
    m_world = NULL;

	//OutputTotal time withouth wating time
	char timeStr[10];
	sprintf(timeStr,"%4.2f\n",m_runSeconds);
	m_timeOutputFile << "Total run time: \t" << timeStr;
	m_timeOutputFile << "---------END OF PROGRAM----------\n";
	m_timeOutputFile.close();
}

void App::message(const String& msg) const {
    renderDevice->clear();
    renderDevice->push2D();
        debugFont->draw2D(renderDevice, msg, renderDevice->viewport().center(), 12, 
            Color3::white(), Color4::clear(), GFont::XALIGN_CENTER, GFont::YALIGN_CENTER);
    renderDevice->pop2D();

    // Force update so that we can see the message
    renderDevice->swapBuffers();
}

GuiTabPane* App::makeFeaturePane(GuiPane* p, featureStat stat, int controlWidth ) {
    debugAssert(notNull(p));

    m_featureTabPane = p->addTabPane();
    for (int j = 0; j < CBFilter::FT_SIZE; j++) {
		String typeStr = (stat) ? String(" variance") : String(" mean");
		m_featureTextures[j][stat] = Texture::fromImage(m_featureNames[j]+typeStr,Image3::createEmpty(m_imgWidth, m_imgHeight));
		p = m_featureTabPane->addTab(m_featureNames[j]);
		m_textureBox[j][stat] = p->addTextureBox(m_featureTextures[j][stat]);
		m_textureBox[j][stat]->setSizeFromInterior(Vector2(float(controlWidth * m_imgWidth / m_imgHeight), float(controlWidth * m_imgWidth / m_imgHeight)));
		m_textureBox[j][stat]->zoomToFit();
		p->pack();
    }
    m_featureTabPane->pack();

    return m_featureTabPane;
}

void App::tick(){
	if ( m_runSeconds == 0)
	{
		///Open log file for timing
		m_timeOutputFile.open("timesLog.txt", std::ios::app);
		m_timeOutputFile << "---------START OF PROGRAM AT " << System::currentTimeString().c_str()  << "---------\n";

		char sceneProperties[50];
		String sceneName;
		switch (m_worldScene)
		{
		case World::CORNELL_BOX:
			sceneName = "Cornell Box";
			break;
		case World::SIBENIK:
			sceneName = "Sibenik Cathedral";
			break;
		case World::CRYTEK_SPONZA:
			sceneName = "Crytek Sponza";
			break;
		}	
		sprintf(sceneProperties,"Scene name: \t%s\n",sceneName.c_str());
		m_timeOutputFile << sceneProperties;
		sprintf(sceneProperties,"Resolution: \t%d x %d\n",m_imgWidth,m_imgHeight);
		m_timeOutputFile << sceneProperties;
		sprintf(sceneProperties,"Sample budget: \t%d\n",m_sampleBudget);
		m_timeOutputFile << sceneProperties;
	}
	m_timer.tick();
}

double App::tock(const std::string & operation){
	m_timer.tock();
	RealTime time = m_timer.elapsedTime();
	char timeStr[20];
	sprintf(timeStr,"%4.2f\n",time);
	m_timeOutputFile << operation << "\t" << timeStr;
	m_runSeconds += time;
	return time;
}

void App::debugTime(const std::string &str){
	m_timeOutputFile << str << std::endl;
}