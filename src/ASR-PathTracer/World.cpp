#include "World.h"

World::World(Scene scene) : m_mode(TRACE) {
    
	Any lightSphere,sibenikCath,sponza;

	begin();
	Stopwatch timer;
	switch (scene){

	case SIBENIK:
		sceneName = "Sibenik_Cathedral";
		frame = CFrame::fromXYZYPRDegrees(-17.4f, -11.0f,   2.2f, -85.3f,  -6.8f,   0.0f);
		FOVDeg = 50.1343;
		//focalDist = 15.0;
		focalDist = 0;
		//blurRadius = 0.15;
		blurRadius = 0;
		ambient = Radiance3(0.0);

		sibenikCath = PARSE_ANY
			( ArticulatedModel::Specification {
				filename = "dabrovic_sibenik/sibenik.zip/sibenik.obj";
				} );

		lightSphere = PARSE_ANY 
			( ArticulatedModel::Specification {
				filename = "sphere.ifs";
				scale = 4.5; 
				preprocess = ( setTwoSided(all(),  false); setMaterial(all(),  UniversalMaterial::Specification { 
					emissive = Color3(15.0,15.0,17.0); }); );
			  });
		insert(ArticulatedModel::create(sibenikCath), Vector3(0, 0, 0));
		insert(ArticulatedModel::create(lightSphere),	CFrame::fromXYZYPRDegrees(6,8.5,0));

		break;

	case CRYTEK_SPONZA:
		sceneName = "Crytek_Sponza";
		frame = CFrame::fromXYZYPRDegrees(  1.1f,   2.2f,   0.0f,  85.2f, -12.6f,   0.0f);
		FOVDeg = 70.1343;
		focalDist = 0.0;
		blurRadius = 0.0;
		ambient = Radiance3::fromARGB(0x4277a6) * 0.0f;

		sponza = PARSE_ANY
			( ArticulatedModel::Specification {
				filename = "crytek_sponza/sponza.zip/sponza.obj";
				scale = 0.005;
			  } );
		lightSphere = PARSE_ANY 
			( ArticulatedModel::Specification {
				filename = "sphere.ifs";
				scale = 7.0; 
				preprocess = ( setTwoSided(all(),  false); setMaterial(all(),  UniversalMaterial::Specification { 
					emissive = Color3(30.0,20.0,10.0); }); ); 
			  });
		insert(ArticulatedModel::create(sponza), Vector3(0, 0, 0));
		insert(ArticulatedModel::create(lightSphere),	CFrame::fromXYZYPRDegrees(12,15,1));
		
		break;

	case CORNELL_BOX:
		sceneName = "Cornell_Box";
		frame = CFrame::fromXYZYPRDegrees( 50.0f,  52.0f,  295.6f,   0.0f,  -2.4f,   0.0f);
		FOVDeg = 29.3811;
		focalDist = 0.0;
		blurRadius = 0.0;	
		ambient = Radiance3(0.0);

		Any leftWall = PARSE_ANY 
			( ArticulatedModel::Specification {
				filename = "square.ifs";
				scale = 200000;
				preprocess = ( setTwoSided(all(),  false); setMaterial(all(),  UniversalMaterial::Specification { 
					lambertian = Color3(0.75,0.25,0.25); }); );
			   });
		Any rightWall = PARSE_ANY 
			( ArticulatedModel::Specification {
				filename = "square.ifs";
				scale = 200000;
				preprocess = ( setTwoSided(all(),  false); setMaterial(all(),  UniversalMaterial::Specification { 
					lambertian = Color3(0.25,0.25,0.75); }); );
			   });
		Any backWall = PARSE_ANY 
			( ArticulatedModel::Specification {
				filename = "square.ifs";
				scale = 200000;
				preprocess = ( setTwoSided(all(),  false); setMaterial(all(),  UniversalMaterial::Specification { 
					lambertian = Color3(0.75,0.75,0.75); }); );
			   });
		Any bottomWall = PARSE_ANY 
			( ArticulatedModel::Specification {
				filename = "square.ifs";
				scale = 200000;
				preprocess = ( setTwoSided(all(),  false); setMaterial(all(),  UniversalMaterial::Specification { 
					lambertian = Color3(0.75,0.75,0.75); }); );
			   });
		Any topWall = PARSE_ANY 
			( ArticulatedModel::Specification {
				filename = "square.ifs";
				scale = 200000;
				preprocess = ( setTwoSided(all(),  false); setMaterial(all(),  UniversalMaterial::Specification { 
					lambertian = Color3(0.75,0.75,0.75); }); );
			   });
		Any mirrorSphere = PARSE_ANY 
			( ArticulatedModel::Specification {
				filename = "sphere.ifs";
				scale = 16.5;
				preprocess =
					( setTwoSided(all(),  false);
						setMaterial(all(), 
								  UniversalMaterial::Specification {
									  //glossy = Color4(0.1f, 0.1f, 0.1f, mirror());
									  lambertian = Color3(.99,.99,.99);
									  //etaTransmit = 1.3f;
									  //etaReflect = 1.0f;
									  //transmissive = Color3(0.2f, 0.5f, 0.7f);
								  });
					  );
			  });
		Any glassSphere = PARSE_ANY 
			( ArticulatedModel::Specification {
				filename = "sphere.ifs";
				scale = 16.5;
				preprocess =
					( setTwoSided(all(),  false);
						setMaterial(all(), 
								  UniversalMaterial::Specification {
									  //glossy = Color4(0.1f, 0.1f, 0.1f, mirror());
									  lambertian = Color3(.99,.99,.99);
									  //etaTransmit = 1.3f;
									  //etaReflect = 1.0f;
									  //transmissive = Color3(0.2f, 0.5f, 0.7f);
								  });
					  );
			  });
		lightSphere = PARSE_ANY 
			( ArticulatedModel::Specification {
				filename = "sphere.ifs";
				scale = 10.0; //CORNELLl Box
				preprocess = ( setTwoSided(all(),  false); setMaterial(all(),  UniversalMaterial::Specification { 
					emissive = Color3(12.0f,12.0f,12.0f); }); ); 
			  });

		insert(ArticulatedModel::create(leftWall),		CFrame::fromXYZYPRDegrees(      1,40.8, 81.6,90));
		insert(ArticulatedModel::create(rightWall),		CFrame::fromXYZYPRDegrees(     99,40.8, 81.6,-90));
		insert(ArticulatedModel::create(backWall),		CFrame::fromXYZYPRDegrees(     50,40.8, 0));
		insert(ArticulatedModel::create(bottomWall),	CFrame::fromXYZYPRDegrees(    50,0.0,81.6,0,-90));
		insert(ArticulatedModel::create(topWall),		CFrame::fromXYZYPRDegrees(    50,81.6,81.6,0,90));
		insert(ArticulatedModel::create(mirrorSphere),	CFrame::fromXYZYPRDegrees(     27,16.5,47));
		insert(ArticulatedModel::create(glassSphere),	CFrame::fromXYZYPRDegrees(     73,16.5,78));
		insert(ArticulatedModel::create(lightSphere),	CFrame::fromXYZYPRDegrees(    50,73,81.6));

		break;
	}
	timer.after("Scene load");

    end();

	for (int L = 0; L < m_surfaceArray.size(); ++L)
		 if (surfaceIsEmissive(m_surfaceArray[L]))
			 lightSurfacesIndex.append(L);
}

inline bool World::surfaceIsEmissive(shared_ptr<Surface> surface) const {
	return ( dynamic_cast<UniversalMaterial* >(surface->resolve().get())->emissive().constant().max() != 0 );
}

Color3 World::getEmissiveFromLighti(int index) const{
	return dynamic_cast<UniversalMaterial* >(m_surfaceArray[lightSurfacesIndex[index] ]->resolve().get())->emissive().constant();
}

void World::getLightSpherei(int index,Sphere& sphere){
	
	sphere.center = m_surfaceArray[lightSurfacesIndex[index] ]->frame().translation;
	AABox box;
	m_surfaceArray[lightSurfacesIndex[index] ]->getObjectSpaceBoundingBox(box);
	sphere.radius = box.extent(1)/2.0;
}


void World::begin() {
    debugAssert(m_mode == TRACE);
    m_surfaceArray.clear();
    m_mode = INSERT;
}


void World::insert(const shared_ptr<ArticulatedModel>& model, const CFrame& frame) {
    Array<shared_ptr<Surface> > posed;
    model->pose(posed, frame);
    for (int i = 0; i < posed.size(); ++i) {
        insert(posed[i]);
    }
}

void World::insert(const shared_ptr<Surface>& m) {
    debugAssert(m_mode == INSERT);
    m_surfaceArray.append(m);
}


void World::end() {
    m_triTree.setContents(m_surfaceArray);
    debugAssert(m_mode == INSERT);
    m_mode = TRACE;
}


bool World::lineOfSight(const Vector3& v0, const Vector3& v1) const {
    debugAssert(m_mode == TRACE);
    
    Vector3 d = v1 - v0;
    float len = d.length();
    Ray ray = Ray::fromOriginAndDirection(v0, d / len);
    float distance = len;
    Tri::Intersector intersector;

    // For shadow rays, try to find intersections as quickly as possible, rather
    // than solving for the first intersection
    static const bool exitOnAnyHit = true, twoSidedTest = true;
    return ! m_triTree.intersectRay(ray, intersector, distance, exitOnAnyHit, twoSidedTest);

}

shared_ptr<Surfel> World::intersect(const Ray& ray, float& distance) const {
    debugAssert(m_mode == TRACE);

    return m_triTree.intersectRay(ray, distance);
}
