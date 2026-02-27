"""
    PatternTransfer{T, I} <: TRANSFER

创建特殊用途的 data Transfer, 用于八叉树层内不同 `rank` 间的辐射函数、配置函数的 Pattern 传输数据.
```
parent::MPIArray{T, IA, 3} where {IA}           待传输的 Pattern `MPIArray`。
reqsIndices::I                                  需求的索引
reqsDatas::Dict{Int, SparseMatrixCSC{T, Int}}   用 SparseMatrixCSC 保存 pole 和 θϕ 方向的数据，Dict建保存 cube 数据。
recv_rk2idcs::Dict{Int, I}                      本 `rank` 接收的 `rank` 以及在本地数据中的索引。
send_rk2idcs::Dict{Int, I}                      本 `rank` 发送的 `rank` 以及在本地数据中的索引。
```
"""
struct PatternTransfer{T, I} <: TRANSFER
    parent::MPIArray{T, IA, 3} where {IA}
    reqsIndices::I
    reqsDatas::Dict{Int, SparseMatrixCSC{T, Int}}
    recv_rk2idcs::Dict{Int, I}
    send_rk2idcs::Dict{Int, I}
end


"""
    PatternTransfer(reqsIndices::NTuple{3, Union{UnitRange{Int}, Vector{Int}}}, a::MPIArray{T, IA, 3}; comm = a.comm, rank = a.myrank, np = MPI.Comm_size(comm)) where {T, IA}

创建缓冲区 `PatternTransfer` 存储辐射函数、配置函数 Pattern 交换数据.
"""
function PatternTransfer(reqsIndices::NTuple{3, Union{UnitRange{Int}, Vector{Int}}}, a::MPIArray{T, IA, 3}; comm = a.comm, rank = a.myrank, np = MPI.Comm_size(comm)) where {T, IA}
    
    # 收集所有进程需求的 indices
    all_reqIndices      =   map(idc -> allgather_VecOrUnitRange(idc; comm = comm), reqsIndices)
    rank2reqIndices     =   Dict(zip(0:(np-1), zip(all_reqIndices...)))
    rank2indices        =   a.rank2indices

    # 需要接收的 rank 和数据的 indices.
    recv_ranks  =   indice2ranks(reqsIndices, rank2indices)
    recv_rank2indices   =   grank2indices(recv_ranks, reqsIndices, rank2indices)

    ## 先创建稀疏矩阵的行列索引并保存在字典里
    # 先找出所有的 cube id
    cubeIndices = Int[]
    for rkcubeIdc in values(recv_rank2indices)
        unique!(sort!(append!(cubeIndices, rkcubeIdc[3])))
    end
    
    cubeIdc2IsJs = Dict([i => (Is = Int[], Js = Int[]) for i in cubeIndices])

    # 其次找出所有 cube 的索引
    for idcs in values(recv_rank2indices)
        for iCube in idcs[3]
            Is, Js = cubeIdc2IsJs[iCube]
            append!(Is, repeat(idcs[1], outer = length(idcs[2])))
            append!(Js, repeat(idcs[2], inner = length(idcs[1])))
        end
    end

    reqsDatas = Dict([i => sparse(  cubeIdc2IsJs[i]..., zeros(T, length(cubeIdc2IsJs[i].Is)), 
                                    size(a, 1), size(a, 2)) for i in cubeIndices])
    
    # 需要发送的 rank 和数据在 a.data 内的 indice
    send_ranks = indice2ranks(a.indices, rank2reqIndices)
    send_rank2indices = remoterank2indices(send_ranks, a.indices, rank2reqIndices)

    PatternTransfer{T, eltype(values(recv_rank2indices))}(a, reqsIndices, reqsDatas, recv_rank2indices, send_rank2indices)

end

"""
    sync!(t::PatternTransfer{T, I}; comm = t.parent.comm, rank = t.parent.myrank, np = MPI.Comm_size(comm)) where{T, I}

同步 `t` 中的数据.
"""
function sync!(t::PatternTransfer{T, I}; comm = t.parent.comm, rank = t.parent.myrank, np = MPI.Comm_size(comm)) where{T, I}

    # parent mpi array
    A = t.parent
    # restoring region
    reqsDatas = t.reqsDatas
    # 缓冲数据，用于接收数据然后再赋值给 reqsDatas 的稀疏矩阵
    buffers =  Dict(k => zeros(T, map(length, v)) for (k, v) in t.recv_rk2idcs)

    # begin sync
    req_all = MPI.Request[]
    begin
        for ghostrank in keys(t.recv_rk2idcs)
            req = MPI.Irecv!(buffers[ghostrank], ghostrank, ghostrank*np + rank, A.comm)
            push!(req_all, req)
        end
        for (remoterank, indices) in t.send_rk2idcs
            req = MPI.Isend(A.data[indices...], remoterank, rank*np + remoterank, A.comm)
            push!(req_all, req)
        end
    end
    MPI.Waitall(MPI.RequestSet(req_all), MPI.Status)
    
    # 将数据存储进
    for (ghostrank, indices) in t.recv_rk2idcs
        buffer = buffers[ghostrank]
        for (i, iCube) in enumerate(indices[3])
            setindex!(reqsDatas[iCube], view(buffer, :, :, i), indices[1], indices[2])
        end
    end
    MPI.Barrier(A.comm)

    return nothing

end
