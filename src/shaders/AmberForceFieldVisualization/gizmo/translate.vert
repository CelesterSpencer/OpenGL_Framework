//============================================================================
// Distributed under the MIT License. Author: Adrian Derstroff
//============================================================================

#version 430

layout(location = 0) in vec4 vertex_position;
out int axisID;

uniform mat4 viewMat;
uniform mat4 projMat;

void main() {
	gl_Position = projMat * viewMat * vertex_position;
	axisID = gl_InstanceID;
}