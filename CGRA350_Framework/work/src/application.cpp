
// std
#include <iostream>
#include <string>
#include <chrono>

// glm
#include <glm/gtc/constants.hpp>
#include <glm/gtc/matrix_transform.hpp>

// project
#include "application.hpp"
#include "cgra/cgra_geometry.hpp"
#include "cgra/cgra_gui.hpp"
#include "cgra/cgra_image.hpp"
#include "cgra/cgra_shader.hpp"
#include "cgra/cgra_wavefront.hpp"

using namespace std;
using namespace cgra;
using namespace glm;

Application::Application(GLFWwindow *window) : m_window(window)
{
    terrain_renderer = make_shared<terrainRenderer>();
    m_cam_pos = vec2(0, 0);
}

void Application::render()
{
    // retrieve the window height
    int width, height;
    glfwGetFramebufferSize(m_window, &width, &height);

    m_windowsize = vec2(width, height); // update window size
    glViewport(0, 0, width, height);    // set the viewport to draw to the entire window
    terrain_renderer->m_windowsize = m_windowsize;

    // clear the back-buffer
    glClearColor(0.3f, 0.3f, 0.4f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // enable flags for normal/forward rendering
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LESS);

    // projection matrix
    mat4 proj = perspective(1.f, float(width) / height, 0.1f, 1000.f);

    // view matrix
    mat4 view = translate(mat4(1), vec3(0, 0, -m_distance)) * rotate(mat4(1), m_pitch, vec3(1, 0, 0)) * rotate(mat4(1), m_yaw, vec3(0, 1, 0));

    // helpful draw options
    //if (m_show_grid)
        //drawGrid(view, proj);
    //if (m_show_axis)
        //drawAxis(view, proj);
    glPolygonMode(GL_FRONT_AND_BACK, (m_showWireframe) ? GL_LINE : GL_FILL);

    // draw
    if (show_terrain)
        terrain_renderer->render(view, proj);
}

void Application::renderGUI()
{   
	
    // setup window
    ImGui::SetNextWindowPos(ImVec2(5, 5), ImGuiSetCond_Once);
    ImGui::SetNextWindowSize(ImVec2(300, 350), ImGuiSetCond_Once);
    ImGui::Begin("Options", 0);

    // display current camera parameters
    ImGui::Text("Application %.3f ms/frame (%.1f FPS)", 1000.0f / ImGui::GetIO().Framerate, ImGui::GetIO().Framerate);
    ImGui::Separator();

    // helpful drawing options
    //ImGui::Checkbox("Show axis", &m_show_axis);
    //ImGui::SameLine();
    //ImGui::Checkbox("Show grid", &m_show_grid);
    //ImGui::Checkbox("Wireframe", &m_showWireframe);
    //ImGui::SameLine();
    //if (ImGui::Button("Screenshot"))
    //    rgba_image::screenshot(true);

    //ImGui::Separator();
    if (ImGui::Checkbox("Terrain", &show_terrain))
    ImGui::Separator();

    if (ImGui::Checkbox("Trees (Canceled)", &show_trees))
        ImGui::Separator();
	
    if (show_terrain){
        if (ImGui::CollapsingHeader("Options")) {
            ImGui::Indent();
            terrain_renderer->renderGUI();
            ImGui::Unindent();
        }
    }

    if (show_trees) {
        
    }

    // finish creating window
    ImGui::End();
}

void Application::cursorPosCallback(double xpos, double ypos){
    vec2 whsize = m_windowsize / 2.0f;

    double y0 = glm::clamp((m_mousePosition.y - whsize.y) / whsize.y, -1.0f, 1.0f);
    double y = glm::clamp((float(ypos) - whsize.y) / whsize.y, -1.0f, 1.0f);
    double dy = -(y - y0);

    double x0 = glm::clamp((m_mousePosition.x - whsize.x) / whsize.x, -1.0f, 1.0f);
    double x = glm::clamp((float(xpos) - whsize.x) / whsize.x, -1.0f, 1.0f);
    double dx = x - x0;
    if (m_leftMouseDown)
    {
        vec2 whsize = m_windowsize / 2.0f;

        // clamp the pitch to [-pi/2, pi/2]
        m_pitch += float(acos(y0) - acos(y));
        m_pitch = float(glm::clamp(m_pitch, -pi<float>() / 2, pi<float>() / 2));

        // wrap the yaw to [-pi, pi]
        m_yaw += float(acos(x0) - acos(x));
        if (m_yaw > pi<float>())
            m_yaw -= float(2 * pi<float>());
        else if (m_yaw < -pi<float>())
            m_yaw += float(2 * pi<float>());
    } else if (m_rightMouseDown) {
        m_distance += dy * 10;
    } else if (m_middleMouseDown) {
        m_cam_pos += vec2(dx, dy) * 10.f;
    }

    // updated mouse position
    m_mousePosition = vec2(xpos, ypos);
}

void Application::mouseButtonCallback(int button, int action, int mods)
{
    (void)mods; // currently un-used

    // capture is left-mouse down
    if (button == GLFW_MOUSE_BUTTON_LEFT)
        m_leftMouseDown = (action == GLFW_PRESS); // only other option is GLFW_RELEASE
}

void Application::scrollCallback(double xoffset, double yoffset) {
    (void)xoffset; // currently un-used
    m_distance *= pow(1.1f, -yoffset);
}


void Application::keyCallback(int key, int scancode, int action, int mods) {
    (void)key, (void)scancode, (void)action, (void)mods; // currently un-used
}


void Application::charCallback(unsigned int c) {
    (void)c; // currently un-used
}
void Application::resize(int width, int height)
{
	
}
