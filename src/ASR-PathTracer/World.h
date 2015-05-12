#ifndef World_h
#define World_h

#include <G3D/G3DAll.h>

class World {
private:

    
    TriTree                         m_triTree;
    CPUVertexArray                  m_cpuVertexArray;
	Array<shared_ptr<Surface> >     m_surfaceArray;
    enum Mode {TRACE, INSERT}       m_mode;

public:

    Array<shared_ptr<Light> >       lightArray;
	Array<int>						lightSurfacesIndex;
	
	
    Color3              ambient;
	// Camera settings
	float				focalDist;
	float				blurRadius;
	CFrame				frame;
	float				FOVDeg;
	String				sceneName;

	enum Scene { CORNELL_BOX, SIBENIK, CRYTEK_SPONZA};


    World(Scene scene);

    /** Returns true if there is an unoccluded line of sight from v0
        to v1.  This is sometimes called the visibilty function in the
        literature.*/
    bool lineOfSight(const Vector3& v0, const Vector3& v1) const;

	/** Returns true if the surface is emissive (has at least one non zero emissive component */
	inline bool surfaceIsEmissive(shared_ptr<Surface> surface) const;

	Color3 getEmissiveFromLighti(int index) const;

	void getLightSpherei(int index,Sphere& sphere);

    void begin();

    void insert(const shared_ptr<ArticulatedModel>& model, const CFrame& frame = CFrame());
    void insert(const shared_ptr<Surface>& m);
    void end();

    /**\brief Trace the ray into the scene and return the first
       surface hit.

       \param ray In world space 

       \param distance On input, the maximum distance to trace to.  On
       output, the distance to the closest surface.

       \return The surfel hit, or NULL if none.
     */
    shared_ptr<Surfel> intersect(const Ray& ray, float& distance) const;
};

#endif
