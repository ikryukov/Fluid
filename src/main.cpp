#include <stdlib.h>
#include <stdio.h>

#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <map>
#include <cmath>

#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <OpenGL/gl.h>
#include <glm/glm.hpp>
#include <glm/ext.hpp>

using namespace glm;

static void error_callback(int error, const char* description)
{
    fputs(description, stderr);
}

static void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
    if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
        glfwSetWindowShouldClose(window, GL_TRUE);
}

static void checkGlError(const char* op) {
    for (GLint error = glGetError(); error; error = glGetError()) {
        //LOGI("after %s() glError (0x%x)\n", op, error);
        fprintf(stderr, "after %s() glError (0x%x)\n", op, error);
    }
}


std::string loadShaderFromFile(std::string &filename)
{
    std::ifstream m_stream(filename.c_str());
    std::string line;
    std::string out = "";
    if (m_stream.fail())
    {
        printf("unable to open file\n");
        return out;
    }
    while(std::getline(m_stream, line, '\n'))
    {
        out += line + '\n';
    }
    m_stream.close();
    return out;
}


GLuint BuildShader(const char* source, GLenum shaderType)
{
    GLuint shaderHandle = glCreateShader(shaderType);
    glShaderSource(shaderHandle, 1, &source, 0);
    glCompileShader(shaderHandle);
    
    GLint compileSuccess;
    glGetShaderiv(shaderHandle, GL_COMPILE_STATUS, &compileSuccess);
    
    if (compileSuccess == GL_FALSE) {
        GLchar messages[256];
        glGetShaderInfoLog(shaderHandle, sizeof(messages), 0, &messages[0]);
        std::cout << messages;
        //LOGI(messages);
        //fprintf(stderr, "%s", messages)
        // TODO
        //exit(1);
    }
    
    return shaderHandle;
}

GLuint BuildProgram(const char* vertexShaderSource, const char* fragmentShaderSource)
{
    GLuint vertexShader = BuildShader(vertexShaderSource, GL_VERTEX_SHADER);
    GLuint fragmentShader = BuildShader(fragmentShaderSource, GL_FRAGMENT_SHADER);
    
    GLuint programHandle = glCreateProgram();
    glAttachShader(programHandle, vertexShader);
    glAttachShader(programHandle, fragmentShader);
    glLinkProgram(programHandle);
    checkGlError("link");
    GLint linkSuccess;
    glGetProgramiv(programHandle, GL_LINK_STATUS, &linkSuccess);
    checkGlError("get program");
    if (linkSuccess == GL_FALSE) {
        GLchar messages[256];
        glGetProgramInfoLog(programHandle, sizeof(messages), 0, &messages[0]);
        std::cout << messages;
        //LOGI(messages);
        //while(1){}
        // TODO
        //exit(1);
    }
    return programHandle;
}

struct Vertex {
    glm::vec3 Position;
    glm::vec3 Normal;
    
    Vertex(glm::vec3 pos, glm::vec3 normal)
    {
        Position = pos;
        Normal = normal;
    }
    
    Vertex() {}
    ~Vertex() {}
};

struct VBO
{
    GLuint m_VAO;
    GLuint m_vertexBuffer;
    GLuint m_verticesCount;
};

class GLRender {
public:
    GLRender() {
        
    }
    
    VBO prepareObjectForRender(std::vector<vec3> points, std::vector<int> indices)
    {
        std::vector<Vertex> vertices(indices.size());
        for (size_t i = 0; i < indices.size();)
        {
            vec3 ab = points[indices[i + 1]] - points[indices[i]];
            vec3 ac = points[indices[i + 2]] - points[indices[i]];
            vec3 normal = cross(ab, ac);
            
            vertices[i + 0].Position= points[indices[i + 0]];
            vertices[i + 0].Normal = normal;
            vertices[i + 1].Position = points[indices[i + 1]];
            vertices[i + 1].Normal = normal;
            vertices[i + 2].Position = points[indices[i + 2]];
            vertices[i + 2].Normal = normal;
            i += 3;
        }
        
        VBO curr;
        glGenVertexArrays(1, &curr.m_VAO);
        glBindVertexArray(curr.m_VAO);
        
        // Create the VBO for the vertices.
        glGenBuffers(1, &curr.m_vertexBuffer);
        glBindBuffer(GL_ARRAY_BUFFER, curr.m_vertexBuffer);
        glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(vertices[0]), &vertices[0], GL_STATIC_DRAW);
        assert(curr.m_vertexBuffer != 0);
        glEnableVertexAttribArray(m_attribPosition);
        glVertexAttribPointer(m_attribPosition, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (GLvoid*) 0);
        glEnableVertexAttribArray(m_attribNormal);
        glVertexAttribPointer(m_attribNormal, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (GLvoid*) offsetof(Vertex, Normal));
        
        glBindVertexArray(0);
        
        curr.m_verticesCount = (GLuint) indices.size();
        
        return curr;
    }
    
    void init(int width, int height)
    {
        glGenRenderbuffers(1, &m_colorRenderbuffer);
        glBindRenderbuffer(GL_RENDERBUFFER, m_colorRenderbuffer);
        glRenderbufferStorage(GL_RENDERBUFFER, GL_RGBA8, width, height);
        checkGlError("color buffer");
        glGenRenderbuffers(1, &m_depthRenderbuffer);
        glBindRenderbuffer(GL_RENDERBUFFER, m_depthRenderbuffer);
        glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT16, width, height);
        checkGlError("depth buffer");
        glGenFramebuffers(1, &m_framebuffer);
        glBindFramebuffer(GL_FRAMEBUFFER, m_framebuffer);
        //	glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, m_colorbufferTexture, 0);
        glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_RENDERBUFFER, m_colorRenderbuffer);
        checkGlError("GL_COLOR_ATTACHMENT0");
        glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, m_depthRenderbuffer);
        checkGlError("GL_DEPTH_ATTACHMENT");
        glBindFramebuffer(GL_FRAMEBUFFER, m_framebuffer);
        glBindFramebuffer(GL_FRAMEBUFFER, 0);
        checkGlError("frame buffer");
        
        glViewport(0, 0, width, height);
        
        glEnable(GL_DEPTH_TEST);
        glEnable(GL_CULL_FACE);
        glCullFace(GL_BACK);
        
        std::string shaderPath = std::string("Simple.vert");
        std::string vertexShaderSource = loadShaderFromFile(shaderPath);
        shaderPath = std::string("Simple.frag");
        std::string fragmentShaderSource = loadShaderFromFile(shaderPath);
        
        m_shaderProgram = BuildProgram(vertexShaderSource.c_str(),
                                       fragmentShaderSource.c_str());
        checkGlError("program");
        
        // set uniforms
        m_uniformProjection = glGetUniformLocation(m_shaderProgram, "Projection");
        checkGlError("get uniform Projection");
        m_uniformView = glGetUniformLocation(m_shaderProgram, "View");
        checkGlError("get uniform View");
        m_uniformTransform = glGetUniformLocation(m_shaderProgram, "Transform");
        checkGlError("get uniform Transform");
        m_uniformColor = glGetUniformLocation(m_shaderProgram, "Color");
        checkGlError("get uniform Color");
        // set attribs
        m_attribPosition = glGetAttribLocation(m_shaderProgram, "Position");
        checkGlError("getAttrib Position");
        m_attribNormal = glGetAttribLocation(m_shaderProgram, "Normal");
        checkGlError("getAttrib Normal");
        
        m_height = height;
        m_with = width;
        m_rotationAngle = 0.0f;
        addBox();
    }
    
    void addBox()
    {
        std::vector<vec3> vertices(8);
        vertices[0] = vec3(1, 1, 1);
        vertices[1] = vec3(-1, 1, 1);
        vertices[2] = vec3(-1, 1, -1);
        vertices[3] = vec3(1, 1, -1);
        
        vertices[4] = vec3(1, -1, -1);
        vertices[5] = vec3(1, -1, 1);
        vertices[6] = vec3(-1, -1, 1);
        vertices[7] = vec3(-1, -1, -1);
        
        int ind[36] =
        {0, 1, 2,
            2, 3, 0,
            0, 3, 4,
            4, 5, 0,
            5, 4, 7,
            7, 6, 5,
            1, 6, 7,
            7, 2, 1,
            6, 1, 0,
            0, 5, 6,
            7, 4, 3,
            7, 3, 2 };
        std::vector<int> indices(ind, ind + 36);
        VBO curr = prepareObjectForRender(vertices, indices);
        m_renderQueue.push_back(curr);
    }
    
    void render()
    {
        renderSolid();
    }
    
    void renderSolid()
    {
        glClearColor(0.5, 0.5, 0.5, 1.0);
        checkGlError("glClearColor");
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        checkGlError("glClear");
        glViewport(0, 0, m_with, m_height);
        checkGlError("glViewport");
        
        glm::mat4 viewMatrix = glm::lookAt(glm::vec3(10, 0, 10), glm::vec3(0, 0, 0), glm::vec3(0, 1, 0));
        glm::mat4 projectionMatrix = glm::perspectiveFov(45.0f, (m_with + 0.0f), m_height + 0.0f, 0.1f, 1000.0f);
        mat4 identity = mat4(1.0f);
        
        glUseProgram(m_shaderProgram);
        checkGlError("glUseProgram(m_shaderProgram);");
        glUniformMatrix4fv(m_uniformProjection, 1, GL_FALSE, glm::value_ptr(projectionMatrix));
        glUniformMatrix4fv(m_uniformView, 1, GL_FALSE, glm::value_ptr(viewMatrix));
        
        vec4 color(0.8, 0.1, 0.4, 1);
        glUniform4fv(m_uniformColor, 1, value_ptr(color));
        
        for (int i = 0; i < m_renderQueue.size(); ++i) {
            glUniformMatrix4fv(m_uniformTransform, 1, GL_FALSE, glm::value_ptr(glm::mat4x4(1.0f)));
            checkGlError("glUniformMatrix4fv");
            glBindVertexArray(m_renderQueue[i].m_VAO);
            glDrawArrays(GL_TRIANGLES, 0, m_renderQueue[i].m_verticesCount);
            checkGlError("glDrawArrays");
        }
    }
    
private:
    int m_with, m_height;
    float m_rotationAngle;
    
    std::vector<VBO> m_renderQueue;
    
    glm::mat4 m_rotation, m_scale, m_translation;
    glm::mat4 m_modelViewMatrix, m_projectionMatrix;
    
    GLuint m_colorRenderbuffer;
    GLuint m_depthRenderbuffer;
    GLuint m_framebuffer;
    
    GLuint m_shaderProgram;
    
    // uniforms
    GLint m_uniformProjection, m_uniformView, m_uniformTransform, m_uniformColor;
    // attribs
    GLint m_attribPosition, m_attribNormal;
};

int main(void)
{
    GLFWwindow* window;
    glfwSetErrorCallback(error_callback);
    if (!glfwInit())
        exit(EXIT_FAILURE);
    glfwDefaultWindowHints();
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);
    
    int width = 1024, height = 768;
    window = glfwCreateWindow(width, height, "Fluid", NULL, NULL);
    if (!window)
    {
        glfwTerminate();
        exit(EXIT_FAILURE);
    }
    // Print the OpenGL version we are using...
    int l_Major = glfwGetWindowAttrib(window, GLFW_CONTEXT_VERSION_MAJOR);
    int l_Minor = glfwGetWindowAttrib(window, GLFW_CONTEXT_VERSION_MINOR);
    printf("OpenGL: %d.%d\n", l_Major, l_Minor);
    
    
    glfwMakeContextCurrent(window);
    glfwSetKeyCallback(window, key_callback);
    
    glewExperimental=true;
    GLenum err = glewInit();
    if (err != GLEW_OK) {
        std::cout << "GLEW init failed: " << glewGetErrorString(err) << "\n";
        exit(1);
    } else {
        std::cout << "Using GLEW " << glewGetString(GLEW_VERSION) << "\n";
    }
    
    int resX = 32, resY = 32, resZ = 32;
    float size = 0.5f/60 * 5;
    
    GLRender render;
    render.init(width, height);
    bool running = true;
    
//    FluidSolver solver(128, 128, 128, 0.5f/60.0f * 5.0f);
    
    while (!glfwWindowShouldClose(window))
    {
        render.render();
        
        glfwSwapBuffers(window);
        
        if (running) 
        {
            //solver.step(0.1);
        }
        
        glfwPollEvents();
    }
    
    glfwDestroyWindow(window);
    glfwTerminate();
    
    exit(EXIT_SUCCESS);
}
