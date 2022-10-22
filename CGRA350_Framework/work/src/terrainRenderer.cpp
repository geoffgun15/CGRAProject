
// std
#include <iostream>
#include <string>
#include <chrono>
#include <algorithm>
#include <random>

// glm
#include <glm/gtc/constants.hpp>
#include <glm/gtc/matrix_transform.hpp>

// project
#include "terrainRenderer.hpp"
#include "cgra/cgra_geometry.hpp"
#include "cgra/cgra_gui.hpp"
#include "cgra/cgra_shader.hpp"
#include "cgra/cgra_wavefront.hpp"


using namespace std;
using namespace terrain;
using namespace glm;


void basic_terrain_model::draw(const glm::mat4& view, const glm::mat4 proj, const vec4 & clip_plane) {
	mat4 modelview = view * modelTransform;

	glUseProgram(shader); // load shader and variables
	glUniformMatrix4fv(glGetUniformLocation(shader, "uProjectionMatrix"), 1, false, value_ptr(proj));
	glUniformMatrix4fv(glGetUniformLocation(shader, "uModelViewMatrix"), 1, false, value_ptr(modelview));
	glUniform3fv(glGetUniformLocation(shader, "uColor"), 1, value_ptr(color));
	glUniform4fv(glGetUniformLocation(shader, "uClipPlane"), 1, value_ptr(clip_plane));
	glUniform1i(glGetUniformLocation(shader, "textureSampler0"), 3);
	glUniform1i(glGetUniformLocation(shader, "textureSampler1"), 4);
	glUniform1i(glGetUniformLocation(shader, "textureSampler2"), 5);
	glUniform1f(glGetUniformLocation(shader, "scale"), scale);
	glUniform1f(glGetUniformLocation(shader, "blendDist"), blendDist);
	glUniform1f(glGetUniformLocation(shader, "transitionHeight1"), transitionHeight1);
	glUniform1f(glGetUniformLocation(shader, "transitionHeight2"), transitionHeight2);

	glActiveTexture(GL_TEXTURE0 + 3);
	glBindTexture(GL_TEXTURE_2D, sandTexture);
	glActiveTexture(GL_TEXTURE0 + 4);
	glBindTexture(GL_TEXTURE_2D, grassTexture);

	mesh.draw(); // draw
}


terrainRenderer::terrainRenderer() {

	cgra::shader_builder sb;
	sb.set_shader(GL_VERTEX_SHADER, CGRA_SRCDIR + std::string("//res//shaders//terrain//color_vert.glsl"));
	sb.set_shader(GL_FRAGMENT_SHADER, CGRA_SRCDIR + std::string("//res//shaders//terrain//color_frag.glsl"));
	GLuint shader = sb.build();

	m_model.shader = shader;
	m_model.color = vec3(0, 1, 0);

	m_model.modelTransform = translate(mat4(1), vec3(-worldSize / 2, 0, -worldSize / 2));
	generateTerrain(numOctaves);
	m_model.scale = scale;
	

	//Bind Sand Texture
	unsigned int sandTexture;
	glGenTextures(1, &sandTexture);
	glActiveTexture(GL_TEXTURE0 + 3);
	glBindTexture(GL_TEXTURE_2D, sandTexture);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	textureImageSand = cgra::rgba_image(CGRA_SRCDIR + std::string("//res//textures//sand_texture.png"));
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, textureImageSand.size.x, textureImageSand.size.y, 0, GL_RGBA, GL_UNSIGNED_BYTE, textureImageSand.data.data());
	m_model.sandTexture = sandTexture;


	//Bind Grass Texture
	unsigned int grassTexture;
	glGenTextures(1, &grassTexture);
	glActiveTexture(GL_TEXTURE0 + 4);
	glBindTexture(GL_TEXTURE_2D, grassTexture);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	textureImageGrass = cgra::rgba_image(CGRA_SRCDIR + std::string("//res//textures//grass_texture.png"));
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, textureImageGrass.size.x, textureImageGrass.size.y, 0, GL_RGBA, GL_UNSIGNED_BYTE, textureImageGrass.data.data());
	m_model.grassTexture = grassTexture;

	//Bind Snow Texture
	unsigned int snowTexture;
	glGenTextures(1, &snowTexture);
	glActiveTexture(GL_TEXTURE0 + 5);
	glBindTexture(GL_TEXTURE_2D, snowTexture);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	textureImagesnow = cgra::rgba_image(CGRA_SRCDIR + std::string("//res//textures//snow_texture.png"));
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, textureImagesnow.size.x, textureImagesnow.size.y, 0, GL_RGBA, GL_UNSIGNED_BYTE, textureImagesnow.data.data());
	m_model.snowTexture = snowTexture;

}

/* Rendering */
void terrainRenderer::render(const glm::mat4& view, const glm::mat4& proj, const vec4& clip_plane) {
	
	//generate mesh
	mesh_builder plane_mb = generatePlane();

	for (int y = 1; y < m_model.heightMap.size() - 1; y++) {
		for (int x = 1; x < m_model.heightMap.size() - 1; x++) {
			int i = (y - 1) * mapSize + (x - 1);

			plane_mb.vertices[i].pos.y = m_model.heightMap[y][x];
			
			//calc normal
			float normX = m_model.heightMap[y][x - 1] / scale - m_model.heightMap[y][x + 1] / scale;
			float normZ = m_model.heightMap[y - 1][x] / scale - m_model.heightMap[y + 1][x] / scale;
			plane_mb.vertices[i].norm = normalize(vec3(normX, 2, normZ));
		}
	}

		
	//generate texture transition offsets
	for (int y = 0; y <= worldSize / squareSize; y++) {
		for (int x = 0; x <= worldSize / squareSize; x++) {
			int i = y * (worldSize / squareSize + 1) + x;
			plane_mb.vertices[i].offset = homogeneousfbm(x * 1, y * 1, 5);
		}
	}

	m_model.mesh = plane_mb.build();

	// draw the model
	m_model.scale = scale;
	m_model.draw(view, proj, clip_plane);
}

/* GUI */
void terrainRenderer::renderGUI() {

	//generated a new seed and terrain
	if (ImGui::Button("Randomize Seed")) {
		genPermutations();
		generateTerrain(numOctaves);
        
	}

	//Full Random Button
	if (ImGui::Button("Full Randomize")) {
		genPermutations();
		fractalType = 3;
		generateTerrain(numOctaves);

	}

	//Base Terrain Options
	if (ImGui::CollapsingHeader("Base Terrain")) {
		ImGui::Indent();

		//chose terrain type
		if (ImGui::Combo("Terrain Type", &fractalType, "Basic (really basic) \0With Smoothing (better for taller) \0HMF (lower base frequency and freq multiplier) \0fullRandom (really ugly)", 4)) {
			generateTerrain(numOctaves);
		}

		//terrain options 
		if (ImGui::SliderFloat("Scale", &scale, 1, 35, "%.0f", 1.0f)) {
			generateTerrain(numOctaves);
		}

		if (ImGui::SliderFloat("Base Frequency", &baseFrequency, 0, 0.2, "%.3f")) {
			generateTerrain(numOctaves);
		}

		if (ImGui::SliderFloat("Frequency Multiplier", &frequencyMultiplier, 1, 3, "%.1f")) {
			generateTerrain(numOctaves);
		}

		if (ImGui::SliderFloat("Amptitude Multiplier", &amtitudeMultiplier, 0, 0.5, "%.2f")) {
			generateTerrain(numOctaves);
		}

		if (ImGui::SliderInt("Num Octaves", &numOctaves, 0, 10)) {
			generateTerrain(numOctaves);
		}

		//extra options for hybrid multifractal terrain type
		if (fractalType == 2) {
			if (ImGui::SliderFloat("Offset", &offset, -1, 1, "%.2f")) {
				generateTerrain(numOctaves);
			}

			if (ImGui::SliderFloat("H", &H, 0, 1, "%.2f")) {
				generateTerrain(numOctaves);
			}
		}

		ImGui::Unindent();

	}

	//Texture Options
	if (ImGui::CollapsingHeader("Texture")) {
		ImGui::Indent();

		ImGui::SliderFloat("Blend Smoothness", &m_model.blendDist, 0, 4, "%.2f");
		ImGui::SliderFloat("Sand Level", &m_model.transitionHeight1, -1, 1, "%.2f");
		ImGui::SliderFloat("Fog Level", &m_model.transitionHeight2, -1, 1, "%.2f");

		ImGui::Unindent();
	}
}


/* Perlin Noise Implementation */

float terrainRenderer::perlinNoise(float x, float y) {

	//get square corner and point in square coords	
	int X = (int)fmod(x, 256); //square coords
	int Y = (int)fmod(y, 256);
	float xf = fmod(x, 1.0f); //point in square coords
	float yf = fmod(y, 1.0f);

	//get hash for each corner
	int TR = permutations[(permutations[(X + 1) % 256] + Y + 1) % 256]; //top right
	int TL = permutations[(permutations[X       % 256] + Y + 1) % 256]; //top left
	int BR = permutations[(permutations[(X + 1) % 256] + Y)     % 256]; //bottom right
	int BL = permutations[(permutations[X       % 256] + Y)     % 256]; //bottom left

	//get const vector for each corner 
	vec2 TR_ConstVec = getCornerVector(TR);
	vec2 TL_ConstVec = getCornerVector(TL);
	vec2 BR_ConstVec = getCornerVector(BR);
	vec2 BL_ConstVec = getCornerVector(BL);

	//get vertor to point for each corner
	vec2 TR_VectorToPoint = vec2(xf - 1.0f, yf - 1.0f);
	vec2 TL_VectorToPoint = vec2(xf,        yf - 1.0f);
	vec2 BR_VectorToPoint = vec2(xf - 1.0f, yf);
	vec2 BL_VectorToPoint = vec2(xf,		yf);

	//get dot product for each corner
	float TR_Val = dot(TR_ConstVec, TR_VectorToPoint);
	float TL_Val = dot(TL_ConstVec, TL_VectorToPoint);
	float BR_Val = dot(BR_ConstVec, BR_VectorToPoint);
	float BL_Val = dot(BL_ConstVec, BL_VectorToPoint);

	//get fade func of point in square
	float u = fade(xf);
	float v = fade(yf);

	//interp to get result
	float interpTop = lerp(u, TL_Val, TR_Val);
	float interpBottom = lerp(u, BL_Val, BR_Val);
	float result = lerp(v, interpBottom, interpTop); //between -1 and 1

	return result;
}

//Gets the constant vector at a corner given its hash from the permutations table which can then be used to get the
//dot product with the vector towards the point in the square. Does this by returning a different vector based on the
//first 4 bits of the hash of the conter.
vec2 terrainRenderer::getCornerVector(int corner) {
	int hash = corner % 4;
	switch (hash) {
	case 0: return vec2(1.0f, 1.0f); break;
	case 1: return vec2(-1.0f, 1.0f); break;
	case 2: return vec2(-1.0f, -1.0f); break;
	case 3: return vec2(1.0f, -1.0f); break;
	default: return vec2(1); break; //should never happen
	}
}

float terrainRenderer::fade(float t) {
	return t * t * t * (t * (t * 6 - 15) + 10);
}

float terrainRenderer::lerp(float x, float p1, float p2) {
	return p1 + x * (p2 - p1);
}

void terrainRenderer::genPermutations() {
	//shuffle the seed (permutation table)
	std::shuffle(permutations, permutations + 256, default_random_engine());
}

/* Terrain Generation */
void terrainRenderer::generateTerrain(int numOctaves) {
	
	//generate height map
	vector<vector<float>> heightMap(mapSize + 2, vector<float>(mapSize + 2, 0));
	for (int y = 0; y < mapSize + 2; y++) {
		for (int x = 0; x < mapSize + 2; x++) {
			if (fractalType == 0) {
				heightMap[y][x] = homogeneousfbm(x * squareSize, y * squareSize, numOctaves) * scale;
			}
			else if (fractalType == 1) {
				heightMap[y][x] = (heterogeneousfbm(x * squareSize, y * squareSize, numOctaves) - 0.5f) * scale;
			}
			else if (fractalType == 2) {
				heightMap[y][x] = (hybridMultifractal(x * squareSize, y * squareSize, numOctaves) - 0.5f) * scale;
			}
			else if (fractalType == 3) {
				heightMap[y][x] = fullrandom(x * squareSize, y * squareSize, numOctaves) * scale;
			}
		}
	}
	m_model.heightMap = heightMap;

	//generate mesh
	mesh_builder plane_mb = generatePlane();

	for (int y = 1; y < heightMap.size()-1; y++) {
		for (int x = 1; x < heightMap[0].size()-1; x++) {
			int i = (y-1) * mapSize + (x-1);

			plane_mb.vertices[i].pos.y = heightMap[y][x];

			//calc normal
			float normX = heightMap[y][x-1]/scale - heightMap[y][x + 1] / scale; //difference in height of previous vertex and next vertex along the x axis
			float normZ = heightMap[y-1][x] / scale - heightMap[y+1][x] / scale; //difference in height of previous vertex and next vertex along the z axis	
			plane_mb.vertices[i].norm = normalize(vec3(normX, 2, normZ));
		}
	}
	

	//generate texture transition offsets
	for (int y = 0; y < mapSize; y++) {
		for (int x = 0; x < mapSize; x++) {
			int i = y*mapSize + x;
			plane_mb.vertices[i].offset = homogeneousfbm(x * 1, y * 1, 5);
		}
	}

	m_model.mesh = plane_mb.build();

}

mesh_builder terrainRenderer::generatePlane() {

	std::vector<vec3> vertices;
	std::vector<vec2> positions;
	std::vector<int> indices;

	for (int y = 0; y < mapSize; y++) {
		for (int x = 0; x < mapSize; x++) {
			//make vertex
			vertices.push_back(vec3(x* squareSize, 0, y* squareSize));
			positions.push_back(vec2(x,y));

			//make triangles (populate index buffer)
			if (x < mapSize-1 && y < mapSize-1) {
				int currentIndex = y * mapSize + x;

				indices.push_back(currentIndex);
				indices.push_back(currentIndex + 1);
				indices.push_back(currentIndex + mapSize);

				indices.push_back(currentIndex + 1);
				indices.push_back(currentIndex + mapSize + 1);
				indices.push_back(currentIndex + mapSize);
			}
		}
	}

	mesh_builder mb;

	for (int i = 0; i < vertices.size(); i++) {
		mb.push_vertex(mesh_vertex{
						vertices[i],
						vec3(0,1,0),
						positions[i],
						0.0f});
	}

	for (int i = 0; i < indices.size(); i++) {
		mb.push_index(indices[i]);
	}

	return mb;
}

/* Terrain Generation Methods */

float terrainRenderer::fullrandom(float x, float y, int numOctaves) {
	float CurrentHeight = 0;

	//add multiple octaves (frequencies that are double the last frequancy and half the amptitude) together to make rough terrain.
	for (int i = 0; i < numOctaves; i++) {
		float r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
		CurrentHeight += r;

	}

	return CurrentHeight;
}


float terrainRenderer::homogeneousfbm(float x, float y, int numOctaves) {
	float CurrentHeight = 0;

	//add multiple octaves (frequencies that are double the last frequancy and half the amptitude) together to make rough terrain.
	for (int i = 0; i < numOctaves; i++) {
		float noise = perlinNoise(x * baseFrequency * pow(frequencyMultiplier, i), y * baseFrequency * pow(frequencyMultiplier, i)) * (float)pow(amtitudeMultiplier, i);

		//add noise to current height function value
		CurrentHeight += noise;

	}

	return CurrentHeight;
}

float terrainRenderer::heterogeneousfbm(float x, float y, int numOctaves) { //FBM with smoothing
	float CurrentHeight = 0;
	float weight = 1.0f;

	//add multiple octaves (frequencies that are double the last frequancy and half the amptitude) together to make rough terrain.
	//weight the amptitude of each frequqncy by the current height of the function to smooth out valleys.
	for (int i = 0; i < numOctaves; i++) {
		float noise = perlinNoise(x*baseFrequency*pow(frequencyMultiplier,i), y*baseFrequency*pow(frequencyMultiplier,i));

		//move range to [0..1] form [-1..1]
		//adding octaves increases the height of mountains instead of just increasing the roughness (average added height is positive instead of 0).
		noise = (noise + 1) / 2.0f;
		
		//reduce amptitude of higher frequencies
		noise *= (float)pow(amtitudeMultiplier, i);
		
		//scale by current height (to smooth valleys)
		float scaledNoise = noise * weight;
		
		//add noise to current height function value
		CurrentHeight += scaledNoise;

		//get the new weight for the next octave
		weight = fmin(1.0f, CurrentHeight);
	}

	return CurrentHeight;
}

float terrainRenderer::hybridMultifractal(float x, float y, int numOctaves) {
	float CurrentHeight = 0;
	float weight = 1.0f;
	float previousNoiseVal = 1.0f;

	//add multiple octaves (frequencies that are double the last frequancy and half the amptitude) together to make rough terrain.
	//weight the amptitude of each frequqncy by the current height of the function to smooth out valleys.
	for (int i = 0; i < numOctaves; i++) {
		float noise = (perlinNoise(x * baseFrequency * pow(frequencyMultiplier, i), y * baseFrequency * pow(frequencyMultiplier, i)) + offset) * pow(pow(amtitudeMultiplier, i), H);

		//scale by previous frequency
		float scaledNoise = noise * weight;

		//add noise to current height function value
		CurrentHeight += scaledNoise;

		//get the new weight for the next octave
		weight = fmin(1.0f, scaledNoise);
	}

	return CurrentHeight;
}














