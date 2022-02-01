#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include <iostream>
#include <string>
#include <list>
#include <vector>
#include <limits>
#include <functional>

#include "Geometry.h"
#include "GLDebug.h"
#include "Log.h"
#include "ShaderProgram.h"
#include "Shader.h"
#include "Texture.h"
#include "Window.h"
#include "Camera.h"

#include "glm/glm.hpp"
#include "glm/gtc/type_ptr.hpp"
#include <glm/gtc/constants.hpp>

bool paused = false;
bool resetPositions = false;


// We gave this code in one of the tutorials, so leaving it here too
void updateGPUGeometry(GPU_Geometry &gpuGeom, CPU_Geometry const &cpuGeom) {
	gpuGeom.bind();
	gpuGeom.setVerts(cpuGeom.verts);
	gpuGeom.setTextureCoordinates(cpuGeom.textureCoordinates);
	gpuGeom.setNormals(cpuGeom.normals);
}

// EXAMPLE CALLBACKS
class Assignment4 : public CallbackInterface {

public:
	Assignment4() : camera(0.0, 0.0, 2.0), aspect(1.0f) {
	}

	virtual void keyCallback(int key, int scancode, int action, int mods) {
		if (key == GLFW_KEY_SPACE && action == GLFW_PRESS)
		{
			paused = !paused;
		}
		if (key == GLFW_KEY_R && action == GLFW_PRESS)
		{
			resetPositions = !resetPositions;
		}
	}
	virtual void mouseButtonCallback(int button, int action, int mods) {
		if (button == GLFW_MOUSE_BUTTON_RIGHT) {
			if (action == GLFW_PRESS) {
				rightMouseDown = true;
			} else if (action == GLFW_RELEASE) {
				rightMouseDown = false;
			}
		}
	}
	virtual void cursorPosCallback(double xpos, double ypos) {
		if (rightMouseDown) {
			double dx = xpos - mouseOldX;
			double dy = ypos - mouseOldY;
			camera.incrementTheta(dy);
			camera.incrementPhi(dx);
		}
		mouseOldX = xpos;
		mouseOldY = ypos;
	}
	virtual void scrollCallback(double xoffset, double yoffset) {
		camera.incrementR(yoffset);
	}
	virtual void windowSizeCallback(int width, int height) {
		// The CallbackInterface::windowSizeCallback will call glViewport for us
		CallbackInterface::windowSizeCallback(width,  height);
		aspect = float(width)/float(height);
	}

	void viewPipeline(ShaderProgram &sp) {
		//rotate the planet to have poles upright
		glm::mat4 rotationFix = glm::mat4(1.0f);
		rotationFix = glm::rotate(rotationFix, glm::radians(90.0f), glm::vec3(1.0, 0.0, 0.0));
		GLint uniMat = glGetUniformLocation(sp, "rotationFix");
		glUniformMatrix4fv(uniMat, 1, GL_FALSE, glm::value_ptr(rotationFix));

		glm::mat4 M = glm::mat4(1.0);
		M = M * rotationFix;//apply the rotation fix
		glm::mat4 V = camera.getView();
		glm::mat4 P = glm::perspective(glm::radians(45.0f), aspect, 0.01f, 1000.f);

		//try rotation fix to the camera position
		glm::vec3 fixedCamPos = rotationFix*glm::vec4(-camera.getPos(),1);
		//std::cout << "updating camPos x: " << camera.getPos().x << std::endl;
		//glUniform3f(glGetUniformLocation(sp, "camPos"), camera.getPos().x, camera.getPos().y, camera.getPos().z);
		glUniform3f(glGetUniformLocation(sp, "camPos"), fixedCamPos.x, fixedCamPos.y, fixedCamPos.z);

		GLint location = glGetUniformLocation(sp, "light");
		//glm::vec3 light = camera.getPos();
		glm::vec3 light = glm::vec3(0.0f, 0.0f, 0.0f);
		//apply the rotate fix to the light position as well too keep consistent with provided locations
		//light = glm::vec4(light,1) * rotationFix;
		//light = rotationFix * glm::vec4(light,1) ;
		//light = glm::vec4(light,1);
		glUniform3fv(location, 1, glm::value_ptr(light));

		uniMat = glGetUniformLocation(sp, "M");
		glUniformMatrix4fv(uniMat, 1, GL_FALSE, glm::value_ptr(M));
		uniMat = glGetUniformLocation(sp, "V");
		glUniformMatrix4fv(uniMat, 1, GL_FALSE, glm::value_ptr(V));
		uniMat = glGetUniformLocation(sp, "P");
		glUniformMatrix4fv(uniMat, 1, GL_FALSE, glm::value_ptr(P));
	}

	glm::mat4 getM() {
		//rotate the planet to have poles upright
		glm::mat4 rotationFix = glm::mat4(1.0f);
		rotationFix = glm::rotate(rotationFix, glm::radians(90.0f), glm::vec3(1.0, 0.0, 0.0));

		glm::mat4 M = glm::mat4(1.0);
		M = M * rotationFix;//apply the rotation fix

		return M;
	}

	glm::mat4 getV() {
		glm::mat4 V = camera.getView();

		return V;
	}

	glm::mat4 getP() {
		glm::mat4 P = glm::perspective(glm::radians(45.0f), aspect, 0.01f, 1000.f);

		return P;
	}

	Camera camera;

private:

	bool rightMouseDown = false;
	float aspect;
	double mouseOldX;
	double mouseOldY;

};

void createPlanet(CPU_Geometry& planet, float planetRadius, float planetResolution, glm::vec3 planetLocation)
{
	//clear geometry, might remove later and clear geometry in the RENDER LOOP
	planet.verts.clear();
	planet.normals.clear();

	for (float theta = 0.0f; theta <= 2*glm::pi<float>(); theta = theta+planetResolution){//step through a 2PI rotation
		for (float phi = 0.0f; phi <= glm::pi<float>(); phi = phi + planetResolution) {
	//for (float theta = 0.0f; theta <= 2*glm::pi<float>(); theta = theta+planetResolution){//step through a 2PI rotation
	//	for (float phi = 0.0f - glm::pi<float>()/2.0f; phi <= glm::pi<float>() - glm::pi<float>()/2.0f; phi = phi + planetResolution) {

			//compute points to be used for each pair of triangles
			glm::vec3 pt1 = glm::vec3(planetRadius * cos(theta) * sin(phi), planetRadius * sin(theta) * sin(phi), planetRadius * cos(phi));
			glm::vec3 pt2 = glm::vec3(planetRadius * cos(theta) * sin(phi+planetResolution), planetRadius * sin(theta) * sin(phi+planetResolution), planetRadius * cos(phi+planetResolution));		
			glm::vec3 pt3 = glm::vec3(planetRadius * cos(theta+planetResolution) * sin(phi), planetRadius * sin(theta+planetResolution) * sin(phi), planetRadius * cos(phi));
			glm::vec3 pt4 = glm::vec3(planetRadius * cos(theta+planetResolution) * sin(phi+planetResolution), planetRadius * sin(theta+planetResolution) * sin(phi+planetResolution), planetRadius * cos(phi+planetResolution));		

			//apply a planetLocation transformation to each vertex to draw the planet in the desired location
			//pt1 = pt1 + planetLocation;
			//pt2 = pt2 + planetLocation;
			//pt3 = pt3 + planetLocation;
			//pt4 = pt4 + planetLocation;

			//place the triangle points into the planet's CPU_Geometry vector
			planet.verts.push_back(pt2);
			planet.verts.push_back(pt1);
			planet.verts.push_back(pt3);
			planet.verts.push_back(pt4);
			planet.verts.push_back(pt2);
			planet.verts.push_back(pt3);

			//compute equivalent texture coordinates for each vertex of the triangles
			glm::vec2 pt1Texture = glm::vec2(theta / (2.0f * glm::pi<float>()), (phi / glm::pi<float>()));
			glm::vec2 pt2Texture = glm::vec2(theta / (2.0f * glm::pi<float>()), ((phi+ planetResolution) / glm::pi<float>()));
			glm::vec2 pt3Texture = glm::vec2((theta + planetResolution) / (2.0f * glm::pi<float>()), (phi / glm::pi<float>()));
			glm::vec2 pt4Texture = glm::vec2((theta + planetResolution) / (2.0f * glm::pi<float>()), ((phi + planetResolution) / glm::pi<float>()));

			//place the equivalent texture coordinates into texture coordinates vertex attributes
			planet.textureCoordinates.push_back(pt2Texture);
			planet.textureCoordinates.push_back(pt1Texture);
			planet.textureCoordinates.push_back(pt3Texture);
			planet.textureCoordinates.push_back(pt4Texture);
			planet.textureCoordinates.push_back(pt2Texture);
			planet.textureCoordinates.push_back(pt3Texture);

			//compute the normal for each vertex of the pair of triangles
			//glm::vec3 pt1Normal = (pt1 - planetLocation) / glm::length(pt1 - planetLocation); //*-1.0f;
			//glm::vec3 pt2Normal = (pt2 - planetLocation) / glm::length(pt2 - planetLocation); //*-1.0f;
			//glm::vec3 pt3Normal = (pt3 - planetLocation) / glm::length(pt3 - planetLocation); //*-1.0f;
			//glm::vec3 pt4Normal = (pt4 - planetLocation) / glm::length(pt4 - planetLocation); //*-1.0f;
			glm::vec3 pt1Normal = (pt1) / glm::length(pt1); //*-1.0f;
			glm::vec3 pt2Normal = (pt2) / glm::length(pt2); //*-1.0f;
			glm::vec3 pt3Normal = (pt3) / glm::length(pt3); //*-1.0f;
			glm::vec3 pt4Normal = (pt4) / glm::length(pt4); //*-1.0f;

			//apply a planetLocation transformation to each vertex to draw the planet normal's in the desired location
			//pt1Normal = pt1Normal + planetLocation;
			//pt2Normal = pt2Normal + planetLocation;
			//pt3Normal = pt3Normal + planetLocation;
			//pt4Normal = pt4Normal + planetLocation;

			//place the planet vertex normals into the planet's CPU_Geometry vector
			planet.normals.push_back(pt2Normal);
			planet.normals.push_back(pt1Normal);
			planet.normals.push_back(pt3Normal);
			planet.normals.push_back(pt4Normal);
			planet.normals.push_back(pt2Normal);
			planet.normals.push_back(pt3Normal);
		}
	}
}

glm::vec3 positionAboutParent(glm::vec3 parentLocation, float distanceToParent, float parentRotationAngle) {


	glm::vec3 positionAboutParent = glm::vec3( distanceToParent*cos(parentRotationAngle), 0.0f,
											  -distanceToParent*sin(parentRotationAngle))+parentLocation;
	return positionAboutParent;

}

glm::vec3 positionAboutParent2(glm::vec3 parentLocation, float distanceToParent) {


	glm::vec3 positionAboutParent = glm::vec3(parentLocation.x+distanceToParent, parentLocation.y, parentLocation.z);
	return positionAboutParent;

}

//this function rotates a planet about it's own axis and incorporates an axial tilt
void rotatePlanet(float& planetRotation, float planetRotRate, ShaderProgram& shader, glm::vec3 planetLocation,
				  float planetAxialTilt, CPU_Geometry& planetCPUGeom) {

	//rotation matrix to rotate the planet about it's polar axis
	glm::mat4 rotM = glm::mat4(1.0f);
	rotM = glm::rotate(rotM, glm::radians(planetRotation), glm::vec3(0.0, 0.0, 1.0));
	//GLint uniMat = glGetUniformLocation(shader, "rotationMatrix");//rotation matrix for relative camera stability in the shader
	//glUniformMatrix4fv(uniMat, 1, GL_FALSE, glm::value_ptr(rotM));

	//transformation matrix to transform the planet to the origin of our scene
	//glm::mat4 trans = glm::mat4(1.0f);
	//trans = glm::translate(trans, -1.0f * planetLocation);

	//transformation matrix to incorporate axil tilt to this planet
	glm::mat4 axialTilt = glm::mat4(1.0f);
	axialTilt = glm::rotate(axialTilt, glm::radians(planetAxialTilt), glm::vec3(0.0, 1.0, 0.0));

	//create a transformation matrix for our planet that incorporates the axial tilt and planet's rotation about its polar axis
	rotM = axialTilt * rotM;
	GLint uniMat = glGetUniformLocation(shader, "rotationMatrix");
	glUniformMatrix4fv(uniMat, 1, GL_FALSE, glm::value_ptr(rotM));

	//create transformation matrix to move the planet to its intended position by the shader
	glm::mat4 trans = glm::mat4(1.0f);
	trans = glm::translate(trans, 1.0f * planetLocation);
	uniMat = glGetUniformLocation(shader, "transformationMatrix");
	glUniformMatrix4fv(uniMat, 1, GL_FALSE, glm::value_ptr(trans));

	//increment the planet rotation progress
	if (paused == false) {

		if (planetRotation + planetRotRate < 360)
			planetRotation = planetRotation + planetRotRate;
		else
			planetRotation = 0.0f;
	}

	//bind to our openGL state machine
	glDrawArrays(GL_TRIANGLES, 0, GLsizei(planetCPUGeom.verts.size()));

}

//this function rotates a planet about it's own axis and incorporates an axial tilt
glm::vec3 animatePlanet(float& planetRotation, float planetRotRate, ShaderProgram& shader, glm::vec3 planetLocation,
	float planetAxialTilt, CPU_Geometry& planetCPUGeom, float& planetOrbit, float planetOrbitalInclination,
	float planetOrbitRate, glm::vec3 parentLocation, std::shared_ptr<Assignment4> a4) {

	//orbit the Earth around the sun
	glm::mat4 orbitTrans = glm::mat4(1.0f);
	orbitTrans = glm::translate(orbitTrans, 1.0f * planetLocation);
	//transformation matrix to incorporate orbital inclination for this planet
	glm::mat4 orbitalInclination = glm::mat4(1.0f);
	orbitalInclination = glm::rotate(orbitalInclination, glm::radians(planetOrbitalInclination), glm::vec3(0.0, 1.0, 0.0));
	glm::mat4 orbitRotM1 = glm::mat4(1.0f);
	glm::mat4 orbitRotM2 = glm::rotate(orbitRotM1, glm::radians(planetOrbit), glm::vec3(0.0, 0.0, 1.0));
	glm::mat4 orbitRotM3 = orbitalInclination*orbitRotM2 * orbitTrans;

	//transformation matrix to fix incorporated  orbital inclination for this planet
	glm::mat4 rotationFix = glm::mat4(1.0f);
	rotationFix = glm::rotate(rotationFix, glm::radians(-planetRotation), glm::vec3(0.0, 0.0, 1.0));
	glm::mat4 axialTiltFix = glm::mat4(1.0f);
	axialTiltFix = glm::rotate(axialTiltFix, glm::radians(-planetAxialTilt), glm::vec3(0.0, 1.0, 0.0));
	glm::mat4 obritalInclinationFix = glm::mat4(1.0f);
	obritalInclinationFix = glm::rotate(obritalInclinationFix, glm::radians(-planetOrbitalInclination), glm::vec3(0.0, 1.0, 0.0));
	//glm::mat4 newPosTrans = orbitRotM;
	//newPlanetLocation = glm::vec3(newPosTrans * glm::vec4(planetLocation, 1));
 
	//rotation matrix to rotate the planet about it's polar axis
	glm::mat4 rotM = glm::mat4(1.0f);
	rotM = glm::rotate(rotM, glm::radians(planetRotation), glm::vec3(0.0, 0.0, 1.0));
	//GLint uniMat = glGetUniformLocation(shader, "rotationMatrix");//rotation matrix for relative camera stability in the shader
	//glUniformMatrix4fv(uniMat, 1, GL_FALSE, glm::value_ptr(rotM));

	//transformation matrix to transform the planet to the origin of our scene
	//glm::mat4 trans = glm::mat4(1.0f);
	//trans = glm::translate(trans, -1.0f * planetLocation);

	//transformation matrix to incorporate axil tilt to this planet
	glm::mat4 axialTilt = glm::mat4(1.0f);
	axialTilt = glm::rotate(axialTilt, glm::radians(planetAxialTilt), glm::vec3(0.0, 1.0, 0.0));


	glm::vec3 axisOfRotation = glm::vec3(glm::rotate(glm::mat4(1.0f), glm::radians(planetAxialTilt),
										glm::vec3(1.0f, 0.0f, 0.0f)) * glm::vec4(0.0f, 1.0f, 0.0f, 0.0f));

	glm::mat4 axialRotationMatrix = glm::rotate(glm::mat4(1.0), glm::radians(planetRotation), axisOfRotation);

	glm::mat4 reverseAxialRotationMatrix = glm::rotate(glm::mat4(1.0), glm::radians(-planetRotation),
											axisOfRotation);

	//create a transformation matrix for our planet that incorporates the axial tilt and planet's rotation about its polar axis
	rotM = orbitRotM3*axialTilt * rotM;
	GLint uniMat = glGetUniformLocation(shader, "rotationMatrix");
	glUniformMatrix4fv(uniMat, 1, GL_FALSE, glm::value_ptr(rotM));

	//rotationFix = axialTiltFix * rotationFix;
	uniMat = glGetUniformLocation(shader, "rotationFix");
	glUniformMatrix4fv(uniMat, 1, GL_FALSE, glm::value_ptr(rotationFix));

	//create transformation matrix to move the planet to its intended position by the shader
	glm::mat4 trans = glm::mat4(1.0f);
	//trans = glm::translate(trans, 1.0f * parentLocation);
	uniMat = glGetUniformLocation(shader, "transformationMatrix");
	glUniformMatrix4fv(uniMat, 1, GL_FALSE, glm::value_ptr(trans));

	if (paused == false) {
		////increment the planet rotation progress
		if (planetRotation + planetRotRate < 360)
			planetRotation = planetRotation + planetRotRate;
		else
			planetRotation = 0.0f;

		//increment the planet orbit progress
		if (planetOrbit + planetOrbitRate < 360)
			planetOrbit = planetOrbit + planetOrbitRate;
		else
			planetOrbit = 0.0f;
	}

	glm::mat4 a4M = a4->getM();
	glm::mat4 a4V = a4->getV();
	glm::mat4 a4P = a4->getP();
	

	glm::vec4 newPlanetLocation = a4P * a4V * a4M * trans * rotM * glm::vec4(planetLocation, 1);

	//bind to our openGL state machine
	glDrawArrays(GL_TRIANGLES, 0, GLsizei(planetCPUGeom.verts.size()));

	//return glm::vec3(rotM * glm::vec4(planetLocation, 1));
	return newPlanetLocation;

}

//this function rotates a planet about it's own axis and incorporates an axial tilt
glm::vec3 animateMoon(float& planetRotation, float planetRotRate, ShaderProgram& shader, glm::vec3 planetLocation,
	float planetAxialTilt, CPU_Geometry& planetCPUGeom, float& planetOrbit, float planetOrbitalInclination,
	float planetOrbitRate, glm::vec3 parentLocation, std::shared_ptr<Assignment4> a4) {

	//orbit the Moon around the earth
	glm::mat4 orbitTrans = glm::mat4(1.0f);
	orbitTrans = glm::translate(orbitTrans, 1.0f * planetLocation);
	//transformation matrix to incorporate orbital inclination for this planet
	glm::mat4 orbitalInclination = glm::mat4(1.0f);
	orbitalInclination = glm::rotate(orbitalInclination, glm::radians(planetOrbitalInclination), glm::vec3(0.0, 1.0, 0.0));
	glm::mat4 orbitRotM1 = glm::mat4(1.0f);
	glm::mat4 orbitRotM2 = glm::rotate(orbitRotM1, glm::radians(planetOrbit), glm::vec3(0.0, 0.0, 1.0));
	glm::mat4 orbitRotM3 = orbitalInclination * orbitRotM2 * orbitTrans;

	//glm::mat4 newPosTrans = orbitRotM;
	//newPlanetLocation = glm::vec3(newPosTrans * glm::vec4(planetLocation, 1));

	//rotation matrix to rotate the planet about it's polar axis
	glm::mat4 rotM = glm::mat4(1.0f);
	rotM = glm::rotate(rotM, glm::radians(planetRotation), glm::vec3(0.0, 0.0, 1.0));
	//GLint uniMat = glGetUniformLocation(shader, "rotationMatrix");//rotation matrix for relative camera stability in the shader
	//glUniformMatrix4fv(uniMat, 1, GL_FALSE, glm::value_ptr(rotM));

	//transformation matrix to transform the planet to the origin of our scene
	//glm::mat4 trans = glm::mat4(1.0f);
	//trans = glm::translate(trans, -1.0f * planetLocation);

	//transformation matrix to incorporate axil tilt to this planet
	glm::mat4 axialTilt = glm::mat4(1.0f);
	axialTilt = glm::rotate(axialTilt, glm::radians(planetAxialTilt), glm::vec3(0.0, 1.0, 0.0));

	//create a transformation matrix for our planet that incorporates the axial tilt and planet's rotation about its polar axis
	rotM = orbitRotM3 * axialTilt * rotM;
	GLint uniMat = glGetUniformLocation(shader, "rotationMatrix");
	glUniformMatrix4fv(uniMat, 1, GL_FALSE, glm::value_ptr(rotM));

	//create transformation matrix to move the planet to its intended position by the shader
	glm::mat4 trans = glm::mat4(1.0f);
	//trans = glm::translate(trans, 1.0f * parentLocation);
	uniMat = glGetUniformLocation(shader, "transformationMatrix");
	glUniformMatrix4fv(uniMat, 1, GL_FALSE, glm::value_ptr(trans));

	//increment the planet rotation progress
	if (planetRotation + planetRotRate < 360)
		planetRotation = planetRotation + planetRotRate;
	else
		planetRotation = 0.0f;

	//increment the planet orbit progress
	if (planetOrbit + planetOrbitRate < 360)
		planetOrbit = planetOrbit + planetOrbitRate;
	else
		planetOrbit = 0.0f;

	glm::mat4 a4M = a4->getM();
	glm::mat4 a4V = a4->getV();
	glm::mat4 a4P = a4->getP();

	
	//glm::vec4 newPlanetLocation = a4P * a4V * a4M * trans * rotM * glm::vec4(planetLocation, 1);
	glm::vec4 newPlanetLocation = glm::vec4(planetLocation, 1);

	//bind to our openGL state machine
	glDrawArrays(GL_TRIANGLES, 0, GLsizei(planetCPUGeom.verts.size()));

	//return glm::vec3(rotM * glm::vec4(planetLocation, 1));
	return newPlanetLocation;

}

void resetSceneObjectLocations(float& sunRotation, glm::vec3 &sunLocation, float& earthRotation,
							   float& earthOrbit, glm::vec3& earthLocation, float& moonRotation,
							   float& moonOrbit, glm::vec3& moonLocation) {
	//reset Sun position
	sunRotation = 0.0f;//changes as the animation progresses

	//reset Earth position
	earthRotation = 0.0f;//changes as the animation progresses
	earthOrbit = 0.0f;//changes as the animation progresses
	earthLocation = positionAboutParent(sunLocation, 0.5, sunRotation);

	//reset Moon Position
	moonRotation = 0.0f;//changes as the animation progresses
	moonOrbit = 0.0f;//chnages as the animation progresses
	moonLocation = positionAboutParent(earthLocation, 0.15, earthRotation);
}

int main() {
	Log::debug("Starting main");

	std::cout << "hello in main" << std::endl;

	// WINDOW
	glfwInit();
	Window window(800, 800, "CPSC 453"); // can set callbacks at construction if desired

	GLDebug::enable();

	// CALLBACKS
	auto a4 = std::make_shared<Assignment4>();
	window.setCallbacks(a4);	

	ShaderProgram shader("shaders/test.vert", "shaders/test.frag");

	//generate the sun sphere
	float sunRadius = 0.15f;
	float sunResolution = 0.1f;
	float sunRotation = 0.0f;//changes as the animation progresses
	float sunRotRate = 0.02f;
	float sunAxialTilt = -10.0f;
	glm::vec3 sunLocation = { 0.0f, 0.0f, 0.0f };
	CPU_Geometry sunCPUGeom;
	GPU_Geometry sunGPUGeom;
	createPlanet(sunCPUGeom, sunRadius, sunResolution, sunLocation);
	updateGPUGeometry(sunGPUGeom, sunCPUGeom);
	Texture sunTexture("textures/2k_sun.jpg", GL_NEAREST);

	//generate the earth sphere
	float earthRadius = 0.06f;
	float earthResolution = 0.1f;
	float earthRotation = 0.0f;//changes as the animation progresses
	float earthOrbit = 0.0f;//changes as the animation progresses
	float earthRotRate = 0.5f;
	float earthOrbitRate = 0.1f;
	float earthAxialTilt = -25.0f;
	float earthOrbitalInclination = 20.0f;
	//glm::vec3 earthLocation = { 0.4f, 0.2f, 0.3f };
	glm::vec3 earthLocation = positionAboutParent(sunLocation, 0.5, sunRotation);
	CPU_Geometry earthCPUGeom;
	GPU_Geometry earthGPUGeom;
	createPlanet(earthCPUGeom, earthRadius, earthResolution, earthLocation);
	updateGPUGeometry(earthGPUGeom, earthCPUGeom);
	Texture earthTexture("textures/2k_earth_daymap.jpg", GL_NEAREST);
	//bind earth to the sun location
	
	//generate the moon sphere
	float moonRadius = 0.04f;
	float moonResolution = 0.1f;
	float moonRotation = 0.0f;//changes as the animation progresses
	float moonOrbit = 0.0f;//chnages as the animation progresses
	float moonRotRate = 0.3f;
	float moonOrbitRate = 0.01f;
	float moonAxialTilt = -15.0f;
	float moonOrbitalInclination = 20.0f;
	//moonOrbit, moonOrbitalInclination, moonOrbitRate
	//glm::vec3 moonLocation = { 0.4f, 0.3f, 0.4f };
	glm::vec3 moonLocation = positionAboutParent(earthLocation, 0.15, earthRotation);
	CPU_Geometry moonCPUGeom;
	GPU_Geometry moonGPUGeom;
	createPlanet(moonCPUGeom, moonRadius, moonResolution, moonLocation);
	updateGPUGeometry(moonGPUGeom, moonCPUGeom);
	Texture moonTexture("textures/2k_moon.jpg", GL_NEAREST);

	//generate the stars sphere
	float starsRadius = 50.0f;
	float starsResolution = 0.1f;
	glm::vec3 starsLocation = { 0.0f, 0.0f, 0.0f };
	CPU_Geometry starsCPUGeom;
	GPU_Geometry starsGPUGeom;
	createPlanet(starsCPUGeom, starsRadius, starsResolution, starsLocation);
	//reverse normal directions so we can see the sphere texture from inside the sphere
	//for (int i = 0; i < starsCPUGeom.normals.size(); i++) {
	//	starsCPUGeom.normals.at(i) = starsCPUGeom.normals.at(i) * -1.0f;
	//}
	updateGPUGeometry(starsGPUGeom, starsCPUGeom);
	Texture starsTexture("textures/2k_stars.jpg", GL_NEAREST);

	// The current CPU_Geometry and GPU_Geometry classes are defined in
	// Geometry.h/Geometry.cpp They will work for this assignment, but for some of
	// the bonuses you may have to modify them.


	// RENDER LOOP
	while (!window.shouldClose()) {

		glfwPollEvents();

		glEnable(GL_LINE_SMOOTH);
		glEnable(GL_FRAMEBUFFER_SRGB);
		glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		glEnable(GL_DEPTH_TEST);
		glDisable(GL_CULL_FACE);

		//glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);		

		shader.use();
		a4->viewPipeline(shader);		

		//reset planet positions
		if (resetPositions == true) {
			//std::cout << "resetting positions" << std::endl;
			resetSceneObjectLocations(sunRotation, sunLocation, earthRotation, earthOrbit, earthLocation, moonRotation, moonOrbit, moonLocation);
			resetPositions = false;
		}

		//draw the planets
		//draw Sun
		sunGPUGeom.bind();
		sunTexture.bind();
		//let the shader know we are drawing Sun or stars
		glUniform1i(glGetUniformLocation(shader, "sunOrStars"), 1);
		//rotate the sun about it's own axis
		rotatePlanet(sunRotation, sunRotRate, shader, sunLocation, sunAxialTilt, sunCPUGeom);

		sunTexture.unbind();

		//draw Earth
		earthGPUGeom.bind();
		earthTexture.bind();
		//let the shader know we are not drawing Sun or stars
		glUniform1i(glGetUniformLocation(shader, "sunOrStars"), 0);
		//rotate the Earth about it's own axis
		//rotatePlanet(earthRotation, earthRotRate, shader, earthLocation, earthAxialTilt, earthCPUGeom);
		//glm::vec3 newEarthLocation = animatePlanet(earthRotation, earthRotRate, shader, earthLocation, earthAxialTilt, earthCPUGeom, earthOrbit, earthOrbitalInclination, earthOrbitRate, sunLocation);
		glm::vec3 newEarthLocation = animatePlanet(earthRotation, earthRotRate, shader, earthLocation, earthAxialTilt, earthCPUGeom, earthOrbit, earthOrbitalInclination, earthOrbitRate, sunLocation, a4);
		//glm::vec3 newEarthLocation = earthLocation;;

		//orbit the Earth around the sun
		//glm::mat4 orbitTrans = glm::mat4(1.0f);
		//orbitTrans = glm::translate(orbitTrans, 1.0f * earthLocation);
		//glm::mat4 orbitRotM = glm::mat4(1.0f);
		//orbitRotM = glm::rotate(orbitRotM, glm::radians(earthOrbit), glm::vec3(0.0, 0.0, 1.0));

		//transformation matrix to incorporate orbital inclination for this planet
		//glm::mat4 orbitalInclination = glm::mat4(1.0f);
		//orbitalInclination = glm::rotate(orbitalInclination, glm::radians(earthOrbitalInclination), glm::vec3(0.0, 1.0, 0.0));

		//create a transformation matrix for our planet that incorporates the axial tilt and planet's rotation about its polar axis
		//orbitRotM = orbitalInclination * orbitRotM;
		//newEarthLocation = glm::vec3(orbitRotM * glm::vec4(earthLocation,1));
		//orbitRotM = orbitRotM * orbitTrans;
		//GLint uniMat = glGetUniformLocation(shader, "rotationMatrix");
		//glUniformMatrix4fv(uniMat, 1, GL_FALSE, glm::value_ptr(orbitRotM));

		//create transformation matrix to move the planet to its intended position by the shader
		//orbitTrans = glm::mat4(1.0f);
		//uniMat = glGetUniformLocation(shader, "transformationMatrix");
		//glUniformMatrix4fv(uniMat, 1, GL_FALSE, glm::value_ptr(orbitTrans));

		//increment the planet rotation progress
		//if (earthOrbit + earthOrbitRate < 360)
		//	earthOrbit = earthOrbit + earthOrbitRate;
		//else
		//	earthOrbit = 0.0f;	
		//bind to our openGL state machine
		//glDrawArrays(GL_TRIANGLES, 0, GLsizei(earthCPUGeom.verts.size()));	
		earthTexture.unbind();


		//draw Moon
		//moonLocation = positionAboutParent2(newEarthLocation, 0.05f);
		//moonLocation = glm::vec3(newEarthLocation.x+0.8, newEarthLocation.y, newEarthLocation.z);
		moonGPUGeom.bind();
		moonTexture.bind();
		//let the shader know we are not drawing Sun or stars
		glUniform1i(glGetUniformLocation(shader, "sunOrStars"), 0);
		//rotate the Moon about it's own axis
		rotatePlanet(moonRotation, moonRotRate, shader, moonLocation, moonAxialTilt, moonCPUGeom);
		//animateMoon(moonRotation, moonRotRate, shader, moonLocation, moonAxialTilt, moonCPUGeom, moonOrbit, moonOrbitalInclination, moonOrbitRate, earthLocation, a4);
		//orbit around the earth

		//manually draw moon as child of earth
		//rotation matrix to rotate the planet about it's polar axis
		glm::mat4 rotM = glm::mat4(1.0f);
		//rotM = glm::rotate(rotM, glm::radians(moonRotation), glm::vec3(0.0, 0.0, 1.0));
		//GLint uniMat = glGetUniformLocation(shader, "rotationMatrix");//rotation matrix for relative camera stability in the shader
		//glUniformMatrix4fv(uniMat, 1, GL_FALSE, glm::value_ptr(rotM));

		//transformation matrix to transform the planet to the origin of our scene
		//glm::mat4 trans = glm::mat4(1.0f);
		//trans = glm::translate(trans, -1.0f * planetLocation);

		//transformation matrix to incorporate axil tilt to this planet
		//glm::mat4 axialTilt = glm::mat4(1.0f);
		//axialTilt = glm::rotate(axialTilt, glm::radians(planetAxialTilt), glm::vec3(0.0, 1.0, 0.0));

		//create a transformation matrix for our planet that incorporates the axial tilt and planet's rotation about its polar axis
		//rotM = axialTilt * rotM;
		//GLint uniMat = glGetUniformLocation(shader, "rotationMatrix");
		//glUniformMatrix4fv(uniMat, 1, GL_FALSE, glm::value_ptr(rotM));

		//create transformation matrix to move the planet to its intended position by the shader
		//glm::mat4 trans = glm::mat4(1.0f);
		//trans = glm::translate(trans, 1.0f * planetLocation);
		//uniMat = glGetUniformLocation(shader, "transformationMatrix");
		//glUniformMatrix4fv(uniMat, 1, GL_FALSE, glm::value_ptr(trans));

		//increment the planet rotation progress
		//if (planetRotation + planetRotRate < 360)
		//	planetRotation = planetRotation + planetRotRate;
		//else
		//	planetRotation = 0.0f;

		//bind to our openGL state machine
		glDrawArrays(GL_TRIANGLES, 0, GLsizei(moonCPUGeom.verts.size()));

		moonTexture.unbind();

		//draw background stars
		starsGPUGeom.bind();
		starsTexture.bind();
		//let the shader know we are drawing Sun or stars
		glUniform1i(glGetUniformLocation(shader, "sunOrStars"), 2);
		//glm::mat4 rotM = glm::mat4(1.0f);
		rotM = glm::mat4(1.0f);
		GLint uniMat = glGetUniformLocation(shader, "rotationMatrix");
		//uniMat = glGetUniformLocation(shader, "rotationMatrix");
		glUniformMatrix4fv(uniMat, 1, GL_FALSE, glm::value_ptr(rotM));
		//uniMat = glGetUniformLocation(shader, "rotationMatrix2");
		//glUniformMatrix4fv(uniMat, 1, GL_FALSE, glm::value_ptr(rotM));
		glm::mat4 trans = glm::mat4(1.0f);
		//trans = glm::mat4(1.0f);
		uniMat = glGetUniformLocation(shader, "transformationMatrix");
		glUniformMatrix4fv(uniMat, 1, GL_FALSE, glm::value_ptr(trans));
		glDrawArrays(GL_TRIANGLES, 0, GLsizei(starsCPUGeom.verts.size()));
		starsTexture.unbind();

		glDisable(GL_FRAMEBUFFER_SRGB); // disable sRGB for things like imgui

		window.swapBuffers();
	}

	glfwTerminate();
	return 0;
}
