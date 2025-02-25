#ifndef RENDER_H
#define RENDER_H

#include <stdio.h>
#include <string.h>
#include <GL/glew.h>
#include <GL/freeglut.h>
#include <iostream>
#include <fstream>
#include <string>

static const char *pVS = R"(
#version 330 core
layout (location = 0) in vec2 pos;
layout (location = 1) in vec2 texcoord;
out vec2 texture_coordinates;
void main() { 
    gl_Position = vec4(pos, 1.0,  1.0);
    texture_coordinates = texcoord; 
})";

static const char *pFS = R"(
#version 330 core
in vec2 texture_coordinates;
out vec4 color;

uniform sampler2D texSampler;
uniform sampler2D texSampler2;

void main() {
    //color  = texture(texSampler, texture_coordinates);
    float value = texture(texSampler, texture_coordinates).r;
    // color = vec4(value, 0, 0, 1);
    color  = texture(texSampler2,  vec2(value, 1));
})";

typedef void (*step_function)(int frame, void *ctx_ptr, float *image_data, double *dtime, bool save);
typedef void (*click_function)(void *ctx_ptr, int, int, int, int);
// Frame Counter
int _FRAME_NUMBER = 0;
int _time = 0, _timebase = 0, _FPS_FRAME = 0, _fps = 0;
double _render_time = 0;

step_function _callback;
click_function _click_callback;
void (*_reset_callback)(void *ctx_ptr);
void *_ctx_ptr;

GLuint shader_program;
GLuint vertexbuffer;
GLuint texbuffer;
GLuint texID;

float points[] = {-1, 1, -1, -1, 1, -1, -1, 1, 1, -1, 1, 1};
float texcoords[] = {0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0};

float *_image_data;
double _time_passed = 0;
int _SIDE;

const int _CTR_LINES = 10;
bool _contour_lines[_CTR_LINES]{};
double _contour_pos[_CTR_LINES]{};
unsigned char _contour_color[_CTR_LINES][3] = {{128, 255, 0}, {0, 255, 0}, {0, 255, 128}, {0, 255, 255}, {0, 128, 255}, {0, 0, 255}, {128, 0, 255}, {255, 0, 255}, {255, 0, 128}, {128, 128, 128}};

const int _PLT_COUNT = 6;
unsigned char *_palette_data[_PLT_COUNT];
int _palette_length[_PLT_COUNT];
GLuint _tex_pal_id[_PLT_COUNT];
int _current_palette = 0;

bool _RUN = false;
bool _SINGLE_RENDER = false;

unsigned char _info_char[128];

void initScene(int side, step_function callback_func, void (*reset_cb)(void *ctx_ptr), void *ctx_ptr,
               click_function click_callback, std::string info) {
    _callback = callback_func;
    _ctx_ptr = ctx_ptr;
    _image_data = new float[side * side];
    _reset_callback = reset_cb;
    _click_callback = click_callback;
    _SIDE = side;

    strcpy((char *)_info_char, info.c_str());

    for (int i = 0; i < _CTR_LINES; i++) {
        _contour_pos[i] = 0.5 / _CTR_LINES * (i + 1);
        // _contour_color[i][0] = 0;
        // // _contour_color[i][1] = 64 + 64.0 / _CTR_LINES * (i + 1);
        // _contour_color[i][1] = 255;
        // _contour_color[i][2] = 255.0 / _CTR_LINES * (i + 1);
    }
}

void resetRenderInfo() {
    glClear(GL_COLOR_BUFFER_BIT);
    _FRAME_NUMBER = 0;
    _time_passed = 0;
    _time = 0, _timebase = 0, _FPS_FRAME = 0, _fps = 0;
    _render_time = 0;
}

// Impostare la palette voluta tenendo conto delle correzzioni
void setPalette(int index) {
    unsigned char values[_palette_length[index] * 3];
    int plen = _palette_length[index];
    std::copy(_palette_data[index], _palette_data[index] + plen * 3, values);
    for (int i = 0; i < _CTR_LINES; i++) {
        if (_contour_lines[i]) {
            values[(int)(_contour_pos[i] * plen) * 3 + 0] = _contour_color[i][0];
            values[(int)(_contour_pos[i] * plen) * 3 + 1] = _contour_color[i][1];
            values[(int)(_contour_pos[i] * plen) * 3 + 2] = _contour_color[i][2];
        }
    }

    glActiveTexture(GL_TEXTURE0 + index + 1);
    glBindTexture(GL_TEXTURE_2D, _tex_pal_id[index]);
    glTextureSubImage2D(_tex_pal_id[index], 0, 0, 0, _palette_length[index], 1, GL_RGB, GL_UNSIGNED_BYTE, values);
}

// Crea lo spazio che gli serve
void loadPaletteFromFile(std::string filepath, int index) {
    std::ifstream fileg(filepath);
    int length;
    fileg >> length;
    _palette_data[index] = new unsigned char[length * 3];
    _palette_length[index] = length;
    for (int i = 0; i < length; i++) {
        int value;
        fileg >> value;
        _palette_data[index][i * 3 + 0] = (unsigned char)value;
        fileg >> value;
        _palette_data[index][i * 3 + 1] = (unsigned char)value;
        fileg >> value;
        _palette_data[index][i * 3 + 2] = (unsigned char)value;
    }
}

void onMouseButton(int button, int state, int x, int y) {
    // // Save the left button state
    if (button == GLUT_LEFT_BUTTON && state == 0) {
        // std::cout << state << "  " << x << "  " << y << std::endl;
        int width = glutGet(GLUT_WINDOW_WIDTH);
        int height = glutGet(GLUT_WINDOW_HEIGHT);

        _click_callback(_ctx_ptr, x, y, width, height);
    }
}

void onNormalKeyboardButton(unsigned char ch, int x, int y) {
    // std::cout << "Premuto " << ch << std::endl;
    if (ch == 'P' || ch == 'p') { //P
        glBindTexture(GL_TEXTURE_2D, texID);
        GLint value;
        glGetTexParameteriv(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, &value);
        // std::cout << "Valore texture: " << value << std::endl;
        if (value == GL_LINEAR) {
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
        } else {
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        }
        glBindTexture(GL_TEXTURE_2D, _tex_pal_id[_current_palette]);
        glutPostRedisplay();

    } else if (ch >= '1' && ch <= '9') { // 1 - 9
        int index = ch - 48 - 1;         // normalizziamo 1-9 in 0-8
        if (index <= _PLT_COUNT) {
            glUseProgram(shader_program);
            GLint sampler2Loc = glGetUniformLocation(shader_program, "texSampler2");
            glUniform1i(sampler2Loc, 1 + index);
            _current_palette = index;
            setPalette(index);
        }
        glutPostRedisplay();
    } else if (ch == 'S' || ch == 's') {
        _RUN = !_RUN;
        glutPostRedisplay();
    } else if (ch == 'R' || ch == 'r') {
        _reset_callback(_ctx_ptr);
        resetRenderInfo();
        glutPostRedisplay();
    } else if (ch == 'N' || ch == 'n') {
        _SINGLE_RENDER = true;
        glutPostRedisplay();
    }
}

void onSpecialKeyboardButton(int key, int x, int y) {
    const int numpads[] = {4, 8, 6, 2, 9, 3, 7, 1, 0, -1, 5};

    // std::cout << "Premuto " << key << std::endl;

    if (key >= 100 && key <= 110 && key != 109) {
        int index = numpads[key - 100];
        if (index < _CTR_LINES) {
            _contour_lines[index] = !_contour_lines[index];
            setPalette(_current_palette);
        }
    }
}

void renderInfo(double time) {
    _time = glutGet(GLUT_ELAPSED_TIME);

    _FPS_FRAME++;
    if (_time - _timebase > 1000) {
        _render_time = ((_time - _timebase) * 1.0) / _FPS_FRAME;
        _fps = 1000 / _render_time;
        _timebase = _time;
        _FPS_FRAME = 0;
        // std::cout << "FPS: " << _fps << "  render time: " << _render_time << std::endl;
    }

    glUseProgram(0);
    glColor3f(0.7f, 0.7f, 0.7f);
    glRasterPos2f(-0.95, -0.97);
    const char *text_fps = "FPS: ";
    glutBitmapString(GLUT_BITMAP_HELVETICA_12, (const unsigned char *)text_fps);

    int fps0 = _fps / 100;
    int fps1 = (_fps - fps0 * 100) / 10;
    int fps2 = _fps - fps0 * 100 - fps1 * 10;
    glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, 48 + fps0);
    glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, 48 + fps1);
    glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, 48 + fps2);
    const char *text_render = "  Render Time: ";
    glutBitmapString(GLUT_BITMAP_HELVETICA_12, (const unsigned char *)text_render);

    int dim = 6;
    char arr[dim]{};
    snprintf(arr, dim, "%.1f", _render_time);
    for (int i = 0; i < dim; i++) {
        if (arr[i] >= 48 && arr[i] <= 48 + 9 || arr[i] == '.') {
            glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, arr[i]);
        } else {
            glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, ' ');
        }
    }
    glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, ' ');
    glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, 'm');
    glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, 's');

    const char *text_time = "    Simulation Time: ";
    glutBitmapString(GLUT_BITMAP_HELVETICA_12, (const unsigned char *)text_time);
    char arr_time[8]{};
    snprintf(arr_time, 8, "%.1f", time);
    for (int i = 0; i < 8; i++) {
        if (arr_time[i] >= 32) {
            glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, arr_time[i]);
        }
    }
    glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, ' ');
    glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, 's');

    glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, ' ');
    glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, ' ');
    glutBitmapString(GLUT_BITMAP_HELVETICA_12, _info_char);
}

static void RenderSceneCB() {

    // glClear(GL_COLOR_BUFFER_BIT);
    glUseProgram(shader_program);

    if (_RUN || _FRAME_NUMBER == 0 || _SINGLE_RENDER) {
        double delta_time = 0;
        _callback(_FRAME_NUMBER, _ctx_ptr, _image_data, &delta_time, false);
        _time_passed += delta_time;
        glutPostRedisplay();

        _SINGLE_RENDER = false;
        _FRAME_NUMBER += 1;
    }

    glTextureSubImage2D(texID, 0, 0, 0, _SIDE, _SIDE, GL_RED, GL_FLOAT, _image_data);
    glDrawArrays(GL_TRIANGLES, 0, 6);

    renderInfo(_time_passed);

    glutSwapBuffers();
}

static void RenderWrapper() {
    // if (_FRAME_NUMBER % 100 == 0) {
    //     std::cout << "Renderizzo il " << _FRAME_NUMBER << "-esimo frame" << std::endl;
    // }
    RenderSceneCB();
}

static void CreateVertexBuffer() {

    glGenTextures(1, &texID);
    glBindTexture(GL_TEXTURE_2D, texID);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_R32F, _SIDE, _SIDE, 0, GL_RED, GL_FLOAT, NULL);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    // glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    // glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

    for (int i = 0; i < _PLT_COUNT; i++) {
        loadPaletteFromFile("./palettes/pal" + std::to_string(i + 1) + ".txt", i);
        glGenTextures(1, &_tex_pal_id[i]);
        glActiveTexture(GL_TEXTURE0 + i + 1);
        glBindTexture(GL_TEXTURE_2D, _tex_pal_id[i]);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB32F, _palette_length[i], 1, 0, GL_RGB, GL_UNSIGNED_BYTE, _palette_data[i]);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    }

    GLint samplerLoc = glGetUniformLocation(shader_program, "texSampler");
    GLint sampler2Loc = glGetUniformLocation(shader_program, "texSampler2");
    glUseProgram(shader_program);
    glUniform1i(samplerLoc, 0);
    glUniform1i(sampler2Loc, 1);

    glActiveTexture(GL_TEXTURE0 + 0);
    glBindTexture(GL_TEXTURE_2D, texID);

    glActiveTexture(GL_TEXTURE0 + 1);
    glBindTexture(GL_TEXTURE_2D, _tex_pal_id[0]);

    glGenBuffers(1, &vertexbuffer);
    glBindBuffer(GL_ARRAY_BUFFER, vertexbuffer);
    glBufferData(GL_ARRAY_BUFFER, sizeof(points), points, GL_STATIC_DRAW);

    glGenBuffers(1, &texbuffer);
    glBindBuffer(GL_ARRAY_BUFFER, texbuffer);
    glBufferData(GL_ARRAY_BUFFER, sizeof(texcoords), texcoords, GL_STATIC_DRAW);

    glEnableVertexAttribArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, vertexbuffer);
    glVertexAttribPointer(0,        // attribute 0
                          2,        // 2 floats per vertex
                          GL_FLOAT, // type
                          GL_FALSE, // normalized?
                          0,        // stride -> tightly packed so make it zero
                          0);       // array buffer offset

    glEnableVertexAttribArray(1);
    glBindBuffer(GL_ARRAY_BUFFER, texbuffer);
    glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 0, 0);
}

static void AddShader(GLuint ShaderProgram, const char *pShaderText, GLenum ShaderType) {
    GLuint ShaderObj = glCreateShader(ShaderType);

    if (ShaderObj == 0) {
        fprintf(stderr, "Error creating shader type %d\n", ShaderType);
        exit(0);
    }

    const GLchar *p[1];
    p[0] = pShaderText;
    GLint Lengths[1];
    Lengths[0] = strlen(pShaderText);
    glShaderSource(ShaderObj, 1, p, Lengths);
    glCompileShader(ShaderObj);
    GLint success;
    glGetShaderiv(ShaderObj, GL_COMPILE_STATUS, &success);
    if (!success) {
        GLchar InfoLog[1024];
        glGetShaderInfoLog(ShaderObj, 1024, NULL, InfoLog);
        fprintf(stderr, "Error compiling shader type %d: '%s'\n", ShaderType, InfoLog);
        exit(1);
    }

    glAttachShader(ShaderProgram, ShaderObj);
}

static void CompileShaders() {
    shader_program = glCreateProgram();

    if (shader_program == 0) {
        fprintf(stderr, "Error creating shader program\n");
        exit(1);
    }

    AddShader(shader_program, pVS, GL_VERTEX_SHADER);
    AddShader(shader_program, pFS, GL_FRAGMENT_SHADER);

    GLint Success = 0;
    GLchar ErrorLog[1024] = {0};

    glLinkProgram(shader_program);
    glGetProgramiv(shader_program, GL_LINK_STATUS, &Success);
    if (Success == 0) {
        glGetProgramInfoLog(shader_program, sizeof(ErrorLog), NULL, ErrorLog);
        fprintf(stderr, "Error linking shader program: '%s'\n", ErrorLog);
        exit(1);
    }

    glValidateProgram(shader_program);
    glGetProgramiv(shader_program, GL_VALIDATE_STATUS, &Success);
    if (!Success) {
        glGetProgramInfoLog(shader_program, sizeof(ErrorLog), NULL, ErrorLog);
        fprintf(stderr, "Invalid shader program: '%s'\n", ErrorLog);
        exit(1);
    }

    glUseProgram(shader_program);
}

void debugMessage(GLenum source, GLenum type, GLuint id, GLenum severity, GLsizei length,
                  const GLchar *message, const void *userParam) {
    std::cout << message << std::endl;
}

void renderScene(int *argcp, char **argv) {
    glutInit(argcp, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA);
    glutInitWindowSize(700, 700);
    glutInitWindowPosition(20, 20);
    glutCreateWindow("Problema 9 - Diffusione del Calore");

    glutDisplayFunc(RenderWrapper);

    // Must be done after glut is initialized!
    GLenum res = glewInit();
    if (res != GLEW_OK) {
        fprintf(stderr, "Error: '%s'\n", glewGetErrorString(res));
        return;
    }
    // printf("vendor: %s\n", glGetString(GL_VENDOR));
    // printf("version: %s\n", glGetString(GL_VERSION));
    // printf("renderer: %s\n", glGetString(GL_RENDERER));

    glEnable(GL_DEBUG_OUTPUT);
    glDebugMessageCallback(debugMessage, NULL);

    glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
    CompileShaders();
    CreateVertexBuffer();
    glutMouseFunc(onMouseButton);
    glutKeyboardFunc(onNormalKeyboardButton);
    glutSpecialFunc(onSpecialKeyboardButton);
    glutMainLoop();
}

#endif