#version 410 core
in vec3 Position;
in vec3 Normal;

uniform mat4 Projection;
uniform mat4 View;
uniform mat4 Transform;
uniform vec4 Color;

out vec4 fColor;
out vec3 fNormal;

void main(void)
{
	fColor = Color;
	fNormal = normalize(Normal);
	gl_Position = Projection * View * Transform * vec4(Position, 1.0);
}
