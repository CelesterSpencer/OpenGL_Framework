//============================================================================
// Distributed under the MIT License. Author: Adrian Derstroff
//============================================================================

#version 430

in int axisID;

void main() {
    gl_FragColor = vec4(1.0, 1.0, 1.0, 1.0);
    if (axisID >= 0)
    {
        gl_FragColor = vec4(1.0, 1.0, 1.0, 1.0);
    }
}
