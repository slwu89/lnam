# class data
using Catlab
using StatsBase

students = [
    "Alice", "Bob", "Carol", "Dan", "Erin", "Eve", "Faith",
    "Frank", "Grace", "Heidi", "Ivan", "Judy", "Mallory",
    "Mike", "Niaj", "Olivia", "Oscar", "Peggy", "Rupert",
    "Sybil", "Trent", "Trudy", "Victor", "Walter", "Wendy"
]

classes = [
    "Economics", "History", "English", "Math", "Biology",
    "Chemistry", "Basket Weaving", "Sports"
]

datamat = zeros(Int, (length(students),length(classes)))

for i in eachindex(students)
    datamat[
        i,
        sample(eachindex(classes), sample(2:4), replace=false)
    ] .= 1
end


@present SchUndirectedNamedBipartiteGraph <: SchUndirectedBipartiteGraph begin
    Name::AttrType
    student::Attr(V₁,Name)
    class::Attr(V₂,Name)
end

@acset_type UndirectedNamedBipartiteGraph(SchUndirectedNamedBipartiteGraph, index=[:src, :tgt]) <: AbstractUndirectedBipartiteGraph

data_2mode = @acset UndirectedNamedBipartiteGraph{String} begin
    V₁ = length(students)
    student = students
    V₂ = length(classes)
    class = classes
    E = sum(datamat)
    src = getindex.(findall(datamat .== 1),1)
    tgt = getindex.(findall(datamat .== 1),2)
end

# to get student "Alice" courses:
data_2mode[incident(data_2mode, "Alice", [:src,:student]),[:tgt,:class]]

to_graphviz(
    data_2mode, 
    node_labels=(:student,:class)
    # graph_attrs=Dict(:dpi=>"72",:size=>"",:ratio=>"fill")
)