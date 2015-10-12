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
#ifdef __APPLE__
	#include <OpenGL/gl.h>
#endif
#include <glm/glm.hpp>
#include <glm/ext.hpp>

using namespace glm;

const int INSTANCE_BUFFER_SIZE = 1 << 10;

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

struct LineVertex
{
    vec3 Position;
};

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

struct RenderItem
{
    int VBO;
    mat4x4 m_transform;
};

class GLRender {
public:
    GLRender() {
        
    }
    
    VBO prepareLineForRender(std::vector<vec3> points, std::vector<int> indices)
    {
        std::vector<LineVertex> vertices(indices.size());
        for (size_t i = 0; i < indices.size();)
        {
            vertices[i + 0].Position= points[indices[i + 0]];
            vertices[i + 1].Position = points[indices[i + 1]];
            i += 2;
        }
        
        VBO curr;
        glGenVertexArrays(1, &curr.m_VAO);
        glBindVertexArray(curr.m_VAO);
        
        // Create the VBO for the vertices.
        glGenBuffers(1, &curr.m_vertexBuffer);
        glBindBuffer(GL_ARRAY_BUFFER, curr.m_vertexBuffer);
        glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(vertices[0]), &vertices[0], GL_STATIC_DRAW);
        assert(curr.m_vertexBuffer != 0);
        glEnableVertexAttribArray(m_attribLinePosition);
        glVertexAttribPointer(m_attribLinePosition, 3, GL_FLOAT, GL_FALSE, sizeof(LineVertex), (GLvoid*) 0);
        
		glBindBuffer(GL_ARRAY_BUFFER, m_buffInstanceMatrixVBO);

        for (int i = 0; i < 4; ++i)
        {
            glEnableVertexAttribArray(m_attribLineInstanceMatrix + i);
            checkGlError("glEnableVertexAttribArray");
            glVertexAttribPointer(m_attribLineInstanceMatrix + i, 4, GL_FLOAT, GL_FALSE, sizeof(mat4x4), (GLvoid*) (0 + i * sizeof(vec4)));
            checkGlError("glVertexAttribPointer");
            glVertexAttribDivisor(m_attribLineInstanceMatrix + i, 1);
            checkGlError("glVertexAttribDivisor");
        }
        
        glBindVertexArray(0);
        
        curr.m_verticesCount = (GLuint) indices.size();
        
        
        return curr;
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
        
        shaderPath = std::string("Line.vert");
        vertexShaderSource = loadShaderFromFile(shaderPath);
        shaderPath = std::string("Line.frag");
        fragmentShaderSource = loadShaderFromFile(shaderPath);
        
        m_lineShaderProgram = BuildProgram(vertexShaderSource.c_str(), fragmentShaderSource.c_str());
        
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
        
        // Line render pass
        m_attribLinePosition = glGetAttribLocation(m_lineShaderProgram, "Position");
        checkGlError("glGetAttribLocation Position");
        
        m_attribLineInstanceMatrix = glGetAttribLocation(m_lineShaderProgram, "InstanceMatrix");
        checkGlError("glGetAttribLocation InstanceMatrix");

        glGenBuffers(1, &m_buffInstanceMatrixVBO);
        checkGlError("glGenBuffers");
        glBindBuffer(GL_ARRAY_BUFFER, m_buffInstanceMatrixVBO);
        checkGlError("glBindBuffer");
        glBufferData(GL_ARRAY_BUFFER, INSTANCE_BUFFER_SIZE * sizeof(mat4x4), NULL, GL_DYNAMIC_DRAW);
        checkGlError("glBufferData");
        
        // set uniforms
        m_uniformLineProjection = glGetUniformLocation(m_lineShaderProgram, "Projection");
        checkGlError("get uniform Line Projection");
        m_uniformLineView = glGetUniformLocation(m_lineShaderProgram, "View");
        checkGlError("get uniform View");
        m_uniformLineColor = glGetUniformLocation(m_lineShaderProgram, "Color");
        checkGlError("get uniform Color");
        
        
        m_height = height;
        m_with = width;
        m_rotationAngle = 0.0f;

        int boxId = addBBox();

		for (int i = 0; i < 10; ++i)
		{
			for (int j = 0; j < 10; ++j)
			{
				for (int k = 0; k < 10; ++k)
				{
					RenderItem ri;
					ri.m_transform = mat4(1.0f);
					ri.m_transform = glm::translate(ri.m_transform, vec3(i * 0.1, j * 0.1, k * 0.1));
					ri.m_transform = glm::scale(ri.m_transform, vec3(0.1, 0.1, 0.1));
					ri.VBO = boxId;
					submitToLineRender(ri);
				}
			}
		}
    }
    
    VBO addBox()
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
        return curr;
    }
    
    int nextId()
    {
        static int id = 0;
        return ++id;
    }
    
    int addBBox()
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
        
        int ind[24] =
        {
            0, 1, 1, 2, 2, 3, 3, 0,
            4, 5, 5, 6, 6, 7, 7, 4,
            0, 5, 1, 6, 2, 7, 3, 4
        };
        std::vector<int> indices(ind, ind + 24);
        VBO curr = prepareLineForRender(vertices, indices);
        int id = nextId();
        m_idToVBO[id] = curr;
        return id;
    }
    
    void render()
    {
        glClearColor(0.5, 0.5, 0.5, 1.0);
        checkGlError("glClearColor");
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        checkGlError("glClear");

//        renderSolid();
        renderLines();
    }
    
    void renderSolid()
    {
        glm::mat4 viewMatrix = glm::lookAt(glm::vec3(0, 0, 4), glm::vec3(0, 0, 0), glm::vec3(0, 1, 0));
        glm::mat4 projectionMatrix = glm::perspectiveFov(45.0f, (m_with + 0.0f), m_height + 0.0f, 0.1f, 1000.0f);
        
        glUseProgram(m_shaderProgram);
        checkGlError("glUseProgram(m_shaderProgram);");
        glUniformMatrix4fv(m_uniformProjection, 1, GL_FALSE, glm::value_ptr(projectionMatrix));
        glUniformMatrix4fv(m_uniformView, 1, GL_FALSE, glm::value_ptr(viewMatrix));
        
    }
    
    void renderLines()
    {
        glm::mat4 viewMatrix = glm::lookAt(glm::vec3(0, 0, 4), glm::vec3(0, 0, 0), glm::vec3(0, 1, 0));
        glm::mat4 projectionMatrix = glm::perspectiveFov(45.0f, (m_with + 0.0f), m_height + 0.0f, 0.1f, 1000.0f);
        mat4 identity = mat4(1.0f);
        
        glUseProgram(m_lineShaderProgram);
        checkGlError("glUseProgram(m_lineShaderProgram);");
		glUniformMatrix4fv(m_uniformLineProjection, 1, GL_FALSE, glm::value_ptr(projectionMatrix));
		glUniformMatrix4fv(m_uniformLineView, 1, GL_FALSE, glm::value_ptr(viewMatrix));
        
        vec4 color(1.0, 1.0, 1.0, 1.0);
        glUniform4fv(m_uniformColor, 1, value_ptr(color));

        
        for (std::map<int, std::vector<RenderItem> >::iterator it = m_renderLineQueue.begin(); it != m_renderLineQueue.end(); ++it)
        {
            std::vector<RenderItem> instances = it->second;
            std::vector<mat4x4> buffInstanceMatrices(instances.size());
            
            VBO currVBO = m_idToVBO[it->first];
            glBindVertexArray(currVBO.m_VAO);
            
            for (int i = 0; i < instances.size(); ++i)
            {
                buffInstanceMatrices[i] = instances[i].m_transform;
            }
            glBindBuffer(GL_ARRAY_BUFFER, m_buffInstanceMatrixVBO);
            glBufferData(GL_ARRAY_BUFFER, instances.size() * sizeof(mat4x4), &buffInstanceMatrices[0], GL_DYNAMIC_DRAW);
            glDrawArraysInstanced(GL_LINES, 0, currVBO.m_verticesCount, (GLsizei) instances.size());
            checkGlError("glDrawArraysInstanced");
            glBindVertexArray(0);
        }
    }
    
    void submitToLineRender(RenderItem& ri)
    {
        if (m_renderLineQueue.count(ri.VBO))
        {
            m_renderLineQueue[ri.VBO].push_back(ri);
        }
        else
        {
            m_renderLineQueue[ri.VBO].push_back(ri);
        
        }
    }
    
private:
    int m_with, m_height;
    float m_rotationAngle;
    
    std::map<int, VBO> m_idToVBO;
    std::vector<RenderItem> m_renderSolidQueue;
    std::map<int, std::vector<RenderItem> > m_renderLineQueue;
    
    glm::mat4 m_rotation, m_scale, m_translation;
    glm::mat4 m_modelViewMatrix, m_projectionMatrix;
    
    GLuint m_colorRenderbuffer;
    GLuint m_depthRenderbuffer;
    GLuint m_framebuffer;
    
    GLuint m_shaderProgram;
    
    GLuint m_lineShaderProgram;
    
    // uniforms
    GLint m_uniformProjection, m_uniformView, m_uniformTransform, m_uniformColor;
    // attribs
    GLint m_attribPosition, m_attribNormal;
    
    // Line render pass attribs
    GLint m_attribLinePosition, m_attribLineInstanceMatrix;
    // Line render pass uniforms
    GLint m_uniformLineProjection, m_uniformLineView, m_uniformLineColor;
    GLuint m_buffInstanceMatrixVBO;
    
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
