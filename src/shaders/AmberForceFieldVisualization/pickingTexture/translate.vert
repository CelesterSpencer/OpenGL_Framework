//============================================================================
// Distributed under the MIT License. Author: Adrian Derstroff
//============================================================================

#version 430

layout(location = 0) in vec4 vertex_position;

out int axesID;

uniform mat4 viewMat;
uniform mat4 projMat;

void main() {
	gl_Position = projMat * viewMat * vertex_position;

	axesID = int(gl_InstanceID);
}