# class data
using Catlab
using DataMigrations
using StatsBase

# sample the student-class data
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

# bipartite (2 mode) schema and data
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

# to get student Alice's courses:
data_2mode[incident(data_2mode, "Alice", [:src,:student]),[:tgt,:class]]

to_graphviz(
    data_2mode, 
    node_labels=(:student,:class)
    # graph_attrs=Dict(:dpi=>"72",:size=>"",:ratio=>"fill")
)

# data migration 1 mode

M = @migration SchLabeledGraph SchUndirectedNamedBipartiteGraph begin
    V => V₁
    E => @join begin
        (s1, s2)::V₁
        (e1, e2)::E
        c::V₂
        src(e1) == s1
        tgt(e1) == c
        src(e2) == s2
        tgt(e2) == c
    end
    # src => e₁ ⋅ src
    # tgt => e₂ ⋅ tgt
end