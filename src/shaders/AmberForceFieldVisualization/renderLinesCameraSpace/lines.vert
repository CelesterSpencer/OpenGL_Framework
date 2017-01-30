//============================================================================
// Distributed under the MIT License. Author: Adrian Derstroff
//============================================================================

#version 430

layout(location = 0) in vec4 vertex_position;

uniform mat4 projMat;

void main() {
	gl_Position = projMat * vertex_position;
}
