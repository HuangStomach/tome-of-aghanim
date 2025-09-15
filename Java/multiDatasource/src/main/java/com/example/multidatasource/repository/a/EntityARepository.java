
package com.example.multidatasource.repository.a;

import com.example.multidatasource.model.a.EntityA;
import org.springframework.data.jpa.repository.JpaRepository;
import org.springframework.stereotype.Repository;

@Repository
public interface EntityARepository extends JpaRepository<EntityA, Long> {
}
